/**
 * Module for Haley-Knott scanone (single QTL mapping by HK regression).
 *
 * NOTE: this module is under development
 */

module qtl.core.scanone_hk;

import std.container;
import qtl.core.primitives;
import std.stdio;
import std.algorithm;
import qtl.core.genotype;
import qtl.plugins.input.read_csv;


/*
 * Performs genome scan using the Haley-Knott regression method
 * (regressing phenotypes on conditional genotype probabilities; the
 * multipoint genotype probabilities have already been calculated in
 * calc.genoprob)
 * 
 * The old R/qtl interface read:
 *
 *   n_ind        Number of individuals
 *   n_pos        Number of marker positions
 *   n_gen        Number of different genotypes
 *   Genoprob     Array of conditional genotype probabilities
 *   Addcov       Matrix of additive covariates: Addcov[cov][ind]
 *   n_addcov     Number of columns of Addcov
 *   Intcov       Number of interactive covariates: Intcov[cov][ind]
 *   n_intcov     Number of columns of Intcov
 *   pheno        Phenotype data, as a vector
 *   nphe         Number of phenotypes
 *   weights      Vector of positive weights, of length n_ind
 *   ind_noqtl    Indicators (0/1) of which individuals should be excluded 
 *                from QTL effects.  
 *   Result       Result matrix of size [n_pos x (nphe)] containing the
 *                LOD scores for each phenotype
 *
 * The new (tentative) interface is:
 *
 * Notes: 
 *
 *   This version skips weights, covariates
 *
 * Check: http://www.netlib.org/lapack/lug/node71.html
 *
 **********************************************************************/

double[] scanone_hk(Ms,Ps,Is,Gs)(in Ms markers, in Ps phenotypes, in Is individuals, in Gs genotypes) 
{
/*
  auto sorted_markers = markers.sorted();
  immutable n_mar = sorted_markers.length;
  immutable n_phe = phenotypes.length;
  immutable n_ind = individuals.length;
  immutable n_gen = genotypes.length;
  immutable n_rss = n_phe;
  immutable n_intcov = 0; // later
  immutable n_addcov = 0; // later
  double[] rss;
  rss.reserve(n_rss);
  double[] tmp_pheno;  
  tmp_pheno.reserve(n_ind * n_phe);
  // design matrix for full model (genotypes + cov)
  immutable ncolx = n_gen + (n_gen-1)*n_intcov+n_addcov;
  // immutable rank = ncolx;
  immutable lwork = 3*ncolx + max(n_ind, n_phe);
  auto dwork = new double[][](200,200);
  // dwork.reserve((2*n_ind+1)*ncolx+lwork+(n_ind+ncolx)*n_phe);
  foreach(ref i; dwork)
    i[] = 0;
  auto singular = dwork;
  auto work = singular[ncolx];
  //auto x = work + lwork;
  //auto x_bk = x + n_ind*ncolx;
  //auto yfit = x_bk + n_ind*ncolx;
  //auto coef = yfit + n_ind*n_phe;

  // for(j=0; j<n_ind; j++) 
  //  for(k=0; k<n_phe; k++)
  //    pheno[j+k*n_ind] *= weights[j];

  
*/
  int  i, j, k, k2, s, rank, info, nrss, lwork, ncolx, ind_idx,
    multivar=0;
  double *dwork, x;
  double *x_bk, singular, yfit, rss, rss_det, work, tmppheno, coef;
  double alpha=1.0;
  auto beta=0.0;
  auto tol=TOL;
  double dtmp;

  /* number of rss's, currently multivar is not used so it's always 0 */
  if( (nphe==1) || (multivar==1) )
    nrss = 1;
  else
    nrss = nphe;

  /* allocate memory */
  rss = new double[nrss];            // rss[nphe]
  tmppheno = new double[n_ind*nphe]; // tmppheno[n_ind][nphe]

  /* number of columns in design matrix X for full model */
  ncolx = n_gen + (n_gen-1)*n_intcov+n_addcov; 
  // 1..n_gen ; 1..n_gen * n_intcov ; n_addcov
  /*ncol0 = n_addcov+1;*/
  rank = ncolx;

  /* allocate space and set things up*/
  /*  x = (double *)R_alloc(n_ind*ncol, sizeof(double));
  coef = (double *)R_alloc(ncol, sizeof(double));
  resid = (double *)R_alloc(n_ind, sizeof(double));
  qty = (double *)R_alloc(n_ind, sizeof(double));
  jpvt = (int *)R_alloc(ncol, sizeof(int));
  qraux = (double *)R_alloc(ncol, sizeof(double));
  work = (double *)R_alloc(2 * ncol, sizeof(double)); */
  lwork = 3*ncolx+ MAX(n_ind, nphe);
  if(multivar == 1)
    dwork = new double[(2*n_ind+1)*ncolx+lwork+(n_ind+nphe+ncolx)*nphe];
    // dwork[2*ind][ncolx (geno+cov)] + [max(ind,phe)][3*ncolx] + [ind+phe+ncollx][phe]
  else
    dwork = new double[(2*n_ind+1)*ncolx+lwork+(n_ind+ncolx)*nphe];

  /* split the memory block */
  singular = dwork;
  work = singular + ncolx;
  x = work + lwork;
  x_bk = x + n_ind*ncolx;
  yfit = x_bk + n_ind*ncolx;
  coef = yfit + n_ind*nphe;
  if(multivar == 1) rss_det = coef + ncolx*nphe;

  /* NULL model is now done in R ********************
     (only do it once!)
  for(j=0; j<n_ind; j++) {
    x[j] = 1.0;
    for(k=0; k<n_addcov; k++) 
      x[j+(k+1)*n_ind] = Addcov[k][j];
  }
  F77_CALL(dqrls)(x, &n_ind, &ncol0, pheno, &ny, &tol, coef, resid,
	          qty, &k, jpvt, qraux, work);
  rss0 = 0.0;
  for(j=0; j<n_ind; j++)  rss0 += (resid[j]*resid[j]);
  Null model is now done in R ********************/

  for(j=0; j<n_ind; j++) 
    for(k=0; k<nphe; k++)
      pheno[j+k*n_ind] *= weights[j];
  /* note: weights are really square-root of weights */

  for(i=0; i<n_pos; i++) { /* loop over positions */
    R_CheckUserInterrupt(); /* check for ^C */

    for(k=0; k<n_ind * ncolx; k++) x[k] = 0.0;

    /* fill up X matrix */
    for(j=0; j<n_ind; j++) {
      if(!ind_noqtl[j]) {
	for(k=0; k<n_gen; k++)
	  x[j+k*n_ind] = Genoprob[k][i][j]*weights[j];
      }
      for(k=0; k<n_addcov; k++)
	x[j+(k+n_gen)*n_ind] = Addcov[k][j]*weights[j];
      if(!ind_noqtl[j]) {
	for(k=0,s=0; k<n_gen-1; k++)
	  for(k2=0; k2<n_intcov; k2++,s++) 
	    x[j+(n_gen+n_addcov+s)*n_ind] = Genoprob[k][i][j]*Intcov[k2][j]*weights[j];
      }
    }

    /* linear regression of phenotype on QTL genotype probabilities */
    /*    F77_CALL(dqrls)(x, &n_ind, &ncol, pheno, &ny, &tol, coef, resid,
		    qty, &k, jpvt, qraux, work);
    */
    /* make a copy of x matrix, we may need it */
    memcpy(x_bk, x, n_ind*ncolx*sizeof(double));
    /* make a copy of phenotypes. I'm doing this because 
       dgelss will destroy the input rhs array */
    memcpy(tmppheno, pheno, n_ind*nphe*sizeof(double));
    /* linear regression of phenotype on QTL genotype probabilities */
    mydgelss (&n_ind, &ncolx, &nphe, x, x_bk, pheno, tmppheno, singular,
      &tol, &rank, work, &lwork, &info);

    /* calculate residual sum of squares */
    if(nphe == 1) {
      /* only one phenotype, this is easier */
      /* if the design matrix is full rank */
      if(rank == ncolx) {
        for (k=rank, rss[0]=0.0; k<n_ind; k++)
          rss[0] += tmppheno[k]*tmppheno[k];
      }
      else {
        /* the desigm matrix is not full rank, this is trouble */
        /* calculate the fitted value */
        matmult(yfit, x_bk, n_ind, ncolx, tmppheno, 1);
        /* calculate rss */
        for (k=0, rss[0]=0.0; k<n_ind; k++)
          rss[0] += (pheno[k]-yfit[k]) * (pheno[k]-yfit[k]);
      }
    }
    else { /* multiple phenotypes */
      if(multivar == 1) {
        /* note that the result tmppheno has dimension n_ind x nphe,
	   the first ncolx rows contains the estimates. */
        for (k=0; k<nphe; k++) 
          memcpy(coef+k*ncolx, tmppheno+k*n_ind, ncolx*sizeof(double));
        /* calculate yfit */
        matmult(yfit, x_bk, n_ind, ncolx, coef, nphe);
        /* calculate residual, put the result in tmppheno */
        for (k=0; k<n_ind*nphe; k++)
          tmppheno[k] = pheno[k] - yfit[k];
        mydgemm(&nphe, &n_ind, &alpha, tmppheno, &beta, rss_det);

        /* calculate the determinant of rss */
        /* do Cholesky factorization on rss_det */
        mydpotrf(&nphe, rss_det, &info);
        for(k=0, rss[0]=1.0;k<nphe; k++)
          rss[0] *= rss_det[k*nphe+k]*rss_det[k*nphe+k];
      } 
      else { /* return rss as a vector */
        if(rank == ncolx) {
          for(k=0; k<nrss; k++) {
            ind_idx = k*n_ind;
            for(j=rank, rss[k]=0.0; j<n_ind; j++) {
              dtmp = tmppheno[ind_idx+j];
              rss[k] += dtmp * dtmp;
            }
          }
        }
        else { /* not full rank, this is troubler */
          /* note that the result tmppheno has dimension n_ind x nphe,
          the first ncolx rows contains the estimates. */
          for (k=0; k<nphe; k++) 
            memcpy(coef+k*ncolx, tmppheno+k*n_ind, ncolx*sizeof(double));
          /* calculate yfit */
          matmult(yfit, x_bk, n_ind, ncolx, coef, nphe);
          /* calculate residual, put the result in tmppheno */
          for (k=0; k<n_ind*nphe; k++)
            tmppheno[k] = pheno[k] - yfit[k];
          /* calculate rss */
          for(k=0; k<nrss; k++) {
            ind_idx = k*n_ind;
            for(j=0, rss[k]=0.0; j<n_ind; j++) {
              dtmp = tmppheno[ind_idx+j];
              rss[k] += dtmp * dtmp;
            }
          }
	  
        }
	
      }
    }
    /* make the result */
    /* log10 likelihood */
    for(k=0; k<nrss; k++) 
       Result[k][i] = (double)n_ind/2.0*log10(rss[k]);

  } /* end loop over positions */
  return null;
}

unittest {
  writeln("Unit test " ~ __FILE__);
}

