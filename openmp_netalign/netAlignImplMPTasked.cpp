
/**
 * @return the best set of matches generated at any step of the code
 * according to the choice of matching
 */
double* netAlignMPTasked(CRS_Mat* S, Vec w, graph* L, Vec li, Vec lj, 
                    double* wperm, netalign_parameters opts, double* objective)
{
    double alpha = opts.alpha;
    double beta = opts.beta;
    double gamma = opts.gamma;
    int iter = opts.maxiter;
    int damping_type = opts.dampingtype;
    
    int nthreads;
    #pragma omp parallel
    {
        nthreads=omp_get_num_threads();
    }
    
    assert(S->nrow()==w.length());

    int size=S->nrow();
    int snz=S->nnz();
    int ns=L->sVertices;
    int nt=L->numVertices-ns;
    
    double damping_mult = 1.;

    double* dt=new double[size];
    double* yt=new double[size];
    double* zt=new double[size];
    double* yt1=new double[size];
    double* zt1=new double[size];
    double* dt1=new double[size];
    double* Fvc=new double[snz];
    double* Stdamp=new double[snz];
    double* p=new double[size];
    double* aw=new double[size];
    double* wi=w.values();
    double* omy=new double[size];
    double* omz=new double[size];
    
    double *stvc=new double[snz];
    double *st1vc = new double[snz];
   
    int saveIter=opts.batchrounding;
    batch_rounding histM(saveIter, li.values(), lj.values(), L, 
        S, w.values(), wperm, alpha, beta, opts.approx ? -1 : 1);

    other_max_computer omax1(li.values(), ns, size, nthreads);
    other_max_computer omax2(lj.values(), nt, size, nthreads);
    
    #pragma omp parallel
    {
        #pragma omp for schedule(static) nowait
        for(int i=0;i<snz;i++) {
            stvc[i]=0.;
            st1vc[i]=0.;
        }
    
        #pragma omp for schedule(static) nowait
        for(int i=0;i<size;i++) {
            aw[i]=alpha*wi[i];
            dt1[i]=0.0;
            yt1[i]=0.0;
            zt1[i]=0.0;
        }
    }
    
    int *ic=S->rowPtr();
    int *perm = build_perm(S);
    // compute split points in the matrix S
    std::vector<int> spSvec(nthreads+1); // split points in the matrix S
    compute_split_points(S, spSvec, nthreads);
    
    int lastiter = 0;
    double timepoints[7] = {0};
    
    assert(damping_type == 2);
    
    for(int t=1;t<=iter;t++)
    {
        lastiter = t;
        
        double time0, time1;

        // Swap data from previous iteration
        time0=timer();

        if(t>1) {
            double *tempptr;
            tempptr = dt1; dt1=dt; dt=tempptr;
            tempptr = yt1; yt1=yt; yt=tempptr;
            tempptr = zt1; zt1=zt; zt=tempptr;
            tempptr = st1vc; st1vc=stvc; stvc = tempptr;
        }
        
        time1 = timer(); timepoints[0] += time1-time0; time0=time1;
        
        // 
        
        #pragma omp parallel num_threads(2)
        {

            // Line 3
            // merged with line 4
            
            //
            // Line 4
            //
            #pragma omp task
            {
                double t0 = timer();
                omp_set_num_threads(nthreads/2);
                
                #ifdef CHUNK
                    #pragma omp parallel for schedule(dynamic,CHUNK)
                #else
                    #pragma omp parallel for schedule(static)
                #endif
                for (int i=0; i<size; ++i) {
                    int h=ic[i]-1;
                    int t=ic[i+1]-1;
                    double sum = 0.; 
                    for(int j=h;j<t;j++) { 
                        // compute the value in F
                        double val = beta+st1vc[perm[j]];
                        if(val<0.) {
                            val=0.;
                        } else {
                            if(val>beta) {
                                val=beta;
                            }
                        } 
                        Fvc[j] = val;
                        sum =+ val;
                    }  
                    dt[i]=sum; // y[t] = sums, z[t] = sums
                }
                timepoints[1] += timer() - t0;
            }
            
            //
            // Line 5, 6
            //
           
            #pragma omp task
            {
                double t0 = timer();
                omp_set_num_threads(nthreads/2);
                if (t>1) {
                    omax2.compute(zt1,omy);
                }
                if(t > 1) {
                    omax1.compute(yt1,omz);
                }
                timepoints[2] += timer() - t0;
            }
            
            // 
            // Precompute dampings
            // 
            #pragma omp task
            {
                double t0 = timer();
                omp_set_num_threads(nthreads/2);
                #pragma omp parallel 
                {
                    #pragma omp for schedule(static) nowait
                    for(int i=0;i<snz;i++) {
                        Stdamp[i]=st1vc[i]+st1vc[perm[i]]-beta;
                    }
                    #pragma omp for schedule(static) nowait
                    for(int i=0;i<size;i++) {
                        p[i]=yt1[i]+zt1[i]-aw[i]+dt1[i];
                    }
                }
                timepoints[5] += timer() - t0;
            }
        }
        
        time1 = timer(); timepoints[3] += time1-time0; time0=time1;
        
        assert(damping_type == 2);
        damping_mult *= gamma;

        //
        // Line 7
        //
        #pragma omp parallel
        {                      
            #ifdef CHUNK
                #pragma omp for schedule(dynamic,CHUNK) nowait
            #else
                #pragma omp for schedule(static)
            #endif
            for(int i=0;i<size;i++) {
                int s=ic[i]-1;
                int t=ic[i+1]-1;
                yt[i] = aw[i] + dt[i] - omy[i] + (1-damping_mult)*p[i];
                zt[i] = aw[i] + dt[i] - omz[i] + (1-damping_mult)*p[i];
                //double val = yt[i]+zt[i]-alpha*w(i)-dt[i];
                double val = aw[i] + dt[i] - omz[i] - omy[i];
                for(int j=s;j<t;j++) {
                    stvc[j]=val-Fvc[j];    
                    // apply damping
                    stvc[j]=stvc[j]+(1.0-damping_mult)*(Stdamp[j]);
                }
            }
        }
        
        time1 = timer(); timepoints[4] += time1-time0; time0=time1; 
                    
        histM.add_heuristic(yt);
        histM.add_heuristic(zt);
        
        time1 = timer(); timepoints[6] += time1-time0; time0=time1;
        
        if (!opts.quiet) {
            cout << "iteration " << t << endl;
        }
    }
    
    double time0 = timer();
    *objective=histM.best_objective();
    double* bestm = histM.best_solution();
    histM.detach_best();
    timepoints[6] += timer()-time0; 
    
    
    if (opts.verbose) {
        cout << "Timing Report: " << endl;
        cout << "     Setup : " << timepoints[0]/lastiter << "s/iter" << endl;
        cout << "    Line 3 : " << timepoints[1]/lastiter << "s/iter" << endl;
        cout << "    Line 4 : " << timepoints[2]/lastiter << "s/iter" << endl;
        cout << "  Line 5/6 : " << timepoints[3]/lastiter << "s/iter" << endl;
        cout << "    Line 7 : " << timepoints[4]/lastiter << "s/iter" << endl;
        cout << "   Damping : " << timepoints[5]/lastiter << "s/iter" << endl;
        cout << "  Rounding : " << timepoints[6]/lastiter << "s/iter" << endl;
    }
    return bestm;    
}
