\name{mvSLOUCH-internal}
\title{Internal mvSLOUCH objects}
\alias{.ape2slouch.mvsl}
\alias{.beginEstimationParams}
\alias{.bm.estim}                                
\alias{.bm.phyl.cov}                             
\alias{.bm.sim.est}                              
\alias{.bm.simulate}                             
\alias{.boot.cleanup}                            
\alias{.boot.extract}                            
\alias{.boot.extract.model.components}           
\alias{.boot.getval}                             
\alias{.calc_species_covariance}                 
\alias{.calc_species_mean}                       
\alias{.calc.cov.bm.mv}                          
\alias{.calc.cov.ouch.mv}                        
\alias{.calc.cov.slouch.mv}                      
\alias{.calc.exptA}                              
\alias{.calc.exptA.2dim}                         
\alias{.calc.integral.evAStevA}                  
\alias{.calc.integral.evAStevA.2dim}             
\alias{.calc.mean.bm.mv}                         
\alias{.calc.mean.ouch.mv}                       
\alias{.calc.mean.slouch.mv}                     
\alias{.calc.phyl.cov}                           
\alias{.calc.phyl.halflife}                      
\alias{.calc.phyl.LogLik.traits}                 
\alias{.calc.vec.dist}                           
\alias{.calcCI}                                  
\alias{.calcESSanalytical}                       
\alias{.calcESSsim}                              
\alias{.calculate_Tree_dists_forCovariance}      
\alias{.CalcVlq}                                 
\alias{.CalcVlq2}                                
\alias{.CalcVlqStat}                             
\alias{.callPCMBase_mvlik}                       
\alias{.changeSigmatoSyy}                        
\alias{.check_is_better}
\alias{.CheckSanityRegimes}                      
\alias{.cleanUpModelParams}                      
\alias{.col.rot}                                 
\alias{.correct_npar0}                           
\alias{.correct_for_diagonalSigns_parscale}
\alias{.correct_phylGLS_response_by_intercept}   
\alias{.correct.names}                           
\alias{.createMeasurementError}                  
\alias{.createRegimes}                           
\alias{.createStartPointsASyyB}                  
\alias{.decompEigenA.S}                          
\alias{.describe.best.model}                     
\alias{.design_matrix_bm}                        
\alias{.design_matrix_construction}              
\alias{.design_matrix_mvslouch}                  
\alias{.design_matrix_ouch}                      
\alias{.detV}                                    
\alias{.do_phylGLSestimation}                    
\alias{.draw.bm}                                 
\alias{.draw.mvslouch}                           
\alias{.draw.ouch}                               
\alias{.drawOneLineage}                          
\alias{.dzeta.matrix.t}                          
\alias{.dzeta.matrix.t1t2}                       
\alias{.dzeta.matrix.t1t2v2}                     
\alias{.dzeta.matrix.tv2}                        
\alias{.dzetaIJ.matrix}                          
\alias{.dzetaKappa.intercept}                    
\alias{.dzetaKappa.matrix}                       
\alias{.dzetaKappa.matrix.simp.t}                
\alias{.dzetaKappa.matrix2}                      
\alias{.dzetaKappa.simp.intercept}               
\alias{.EvaluatePoint}                           
\alias{.extract_from_signs}                     
\alias{.extract_GLS_results}                     
\alias{.extract_GLS_results_bm}                  
\alias{.extract_GLS_results_mvslouch}            
\alias{.extract_GLS_results_ouch}                
\alias{.extract.model.key.properties}            
\alias{.generate.all.model.setups}               
\alias{.generate.basic.model.setups}             
\alias{.generate.ext.model.setups}               
\alias{.generate.fund.model.setups}              
\alias{.generate.list.of.model.setups}           
\alias{.generate.model.description}              
\alias{.generate.model.se}                       
\alias{.generate.univ.model.setups}              
\alias{.generatenames}                           
\alias{.getESScriteria}                          
\alias{.getMVSLphylCovMatrix}                    
\alias{.get_fullGLSintercept}
\alias{.get_fullGLSintercept_mvslouch}
\alias{.getMVSLproc_size}                        
\alias{.givens.mult.M}                           
\alias{.givens.real}                             
\alias{.glsgc.estim}                             
\alias{.InitialRegimeSetup}                      
\alias{.internal_BrownianMotionModel}            
\alias{.internal_mvslouchModel}                  
\alias{.internal_ouchModel}                      
\alias{.inv.rho.givens.real}                     
\alias{.is0_Merror}
\alias{.Kappa.matrix}                            
\alias{.Kappa.matrix.NA}                         
\alias{.M.mult.givens}                           
\alias{.make.states.mvsl}                        
\alias{.maxlik.estim}                            
\alias{.MinusPhylLogLikFunc}                     
\alias{.MinusPhylLogLikFuncMVslouch}             
\alias{.MinusPhylLogLikFuncouch}                 
\alias{.ModelSetupPCMBase}                       
\alias{.mvslouch_to_ouch_model}                  
\alias{.mvslouch.phyl.cov}                       
\alias{.mvslouch.simulate}                       
\alias{.norm.max}                                
\alias{.num.inv.transform.positive}              
\alias{.num.transform.positive}                  
\alias{.oubm.sim.est}                            
\alias{.ouch.phyl.cov}                           
\alias{.ouch.simulate}                           
\alias{.ouch2slouch.mvsl}                        
\alias{.ouou.sim.est}                            
\alias{.OUtoBM_pcmbase_model_box}                
\alias{.packageName}                             
\alias{.par_transform_withPCMBase}               
\alias{.par.inv.transform}                       
\alias{.par.inv.transform.decomp.matrix}         
\alias{.par.inv.transform.decomp.neg.real.matrix}
\alias{.par.inv.transform.decomp.pos.real.matrix}
\alias{.par.inv.transform.decomp.real.matrix}    
\alias{.par.inv.transform.decomp.sym.matrix}     
\alias{.par.inv.transform.invert.matrix}         
\alias{.par.inv.transform.invert.qr.matrix}      
\alias{.par.inv.transform.invert.svd.matrix}     
\alias{.par.inv.transform.lowertri.matrix}       
\alias{.par.inv.transform.orth.matrix.cayley}    
\alias{.par.inv.transform.orth.matrix.givens}    
\alias{.par.inv.transform.symmetric.matrix}      
\alias{.par.inv.transform.twobytwo.matrix}       
\alias{.par.inv.transform.uppertri.matrix}       
\alias{.par.transform}                          
\alias{.par.transform.decomp.matrix}             
\alias{.par.transform.decomp.neg.real.matrix}    
\alias{.par.transform.decomp.pos.real.matrix}    
\alias{.par.transform.decomp.real.matrix}        
\alias{.par.transform.decomp.sym.matrix}         
\alias{.par.transform.invert.matrix}             
\alias{.par.transform.invert.qr.matrix}          
\alias{.par.transform.invert.svd.matrix}         
\alias{.par.transform.lowertri.matrix}           
\alias{.par.transform.orth.matrix.cayley}        
\alias{.par.transform.orth.matrix.givens}        
\alias{.par.transform.symmetric.matrix}          
\alias{.par.transform.twobytwo.matrix}           
\alias{.par.transform.uppertri.matrix}           
\alias{.Params.summary}                          
\alias{.params.summary.bm}                       
\alias{.params.summary.mvslouch}                 
\alias{.params.summary.ouch}                     
\alias{.pcmbaseDphylGaussian_RSS}                
\alias{.pcmbaseDphylOU_GLS}                      
\alias{.PhyloSDEestim}                           
\alias{.phyltree_remove_path_fields}
\alias{.phyltree_remove_tips}
\alias{.postproc_estres}
\alias{.return_best_model}
\alias{.row.rot}
\alias{.set_1var_pcmbase_model_box}
\alias{.set_mean0_pcmbase_model_box}
\alias{.set.estimparams}
\alias{.set.maxiter}
\alias{.set.paramatrizationnames}
\alias{.set_parscale}
\alias{.set.tol}
\alias{.simulVasicekProcPhylTree}
\alias{.SummarizeFullPoint}
\alias{.sym.par}
\alias{.sym.unpar}
\alias{.update_pcmbase_box_params}
\alias{.update_pcmbase_box_params_bm}
\alias{.update_pcmbase_box_params_mvslouch}
\alias{.update_pcmbase_box_params_ouch}
\alias{.internal_phyltree_paths_BL}
\alias{.check_input_trait_data}
\alias{.is_det0}
\alias{.my_cov2cor}
\alias{.my_stop}
\alias{.my_message}
\alias{.my_warning}
\alias{.set_sing_bl_threshold}
\alias{.my_chol}
\alias{.set_pcmbase_model_box_X0}
\alias{.matrixcalc_is.diagonal.matrix} 
\alias{.matrixcalc_is.positive.definite}
\alias{.matrixcalc_is.positive.semi.definite}
\alias{.matrixcalc_is.singular.matrix}
\alias{.matrixcalc_is.square.matrix}
\alias{.matrixcalc_is.symmetric.matrix}
\alias{.calc.phyl.mean}
\alias{.bm.phyl.mean}
\alias{.ouch.phyl.mean}
\alias{.mvslouch.phyl.mean}
\alias{.get_phyl_mean}
\alias{.OUphylreg_internal}
\alias{.id_evolmodel}
\alias{.f_warningsA}
\description{Internal mvSLOUCH objects.}
\details{These are not to be called by the user. They are located in the code files as follows
.ape2slouch.mvsl: trees2slouch.mvsl.R
.beginEstimationParams: modelparams.R
.bm.estim: estimBM.R
.bm.phyl.cov: sdecovariancephyl.R
.bm.sim.est: bootstrap.R
.bm.simulate: simulVasicekproc.R
.boot.cleanup: bootstrap.R
.boot.extract: bootstrap.R
.boot.extract.model.components: bootstrap.R
.boot.getval: bootstrap.R
.calc_species_covariance: sdemoments.R
.calc_species_mean: sdemoments.R
.calc.cov.bm.mv: sdemoments.R
.calc.cov.ouch.mv: sdemoments.R
.calc.cov.slouch.mv: sdemoments.R
.calc.exptA: matrixexps.R
.calc.exptA.2dim: matrixexps.R
.calc.integral.evAStevA: matrixexps.R
.calc.integral.evAStevA.2dim: matrixexps.R
.calc.mean.bm.mv: sdemoments.R
.calc.mean.ouch.mv: sdemoments.R
.calc.mean.slouch.mv: sdemoments.R
.calc.phyl.cov: sdecovariancephyl.R
.calc.phyl.halflife: modelparamssummary.R
.calc.phyl.LogLik.traits: loglik.R
.calc.vec.dist: phylgls.R
.calcCI: ci.R
.calcESSanalytical: getESS.R
.calcESSsim: getESS.R
.calculate_Tree_dists_forCovariance: getESS.R
.CalcVlq: precalcs.R
.CalcVlq2: precalcs.R
.CalcVlqStat: precalcs.R
.callPCMBase_mvlik: loglik.R
.changeSigmatoSyy: evolmodelest.R
.CheckSanityRegimes: regimes.R
.check_is_better: evolmodelest.R
.cleanUpModelParams: modelparamssummary.R
.col.rot: matrixparametrizations.R
.correct_npar0: modelparamssummary.R
.correct_for_diagonalSigns_parscale: PhyloSDEestim.R
.correct_phylGLS_response_by_intercept: phylgls.R
.correct.names: PhyloSDEestim.R
.createMeasurementError: PhyloSDEestim.R
.createRegimes: regimes.R
.createStartPointsASyyB: evolmodelest.R
.decompEigenA.S: precalcs.R
.describe.best.model: evolmodelest.R
.design_matrix_bm: phylgls.R
.design_matrix_construction: phylgls.R
.design_matrix_mvslouch: phylgls.R
.design_matrix_ouch: phylgls.R
.detV: phylgls.R
.do_phylGLSestimation: phylgls.R
.draw.bm: simulVasicekproc.R
.draw.mvslouch: simulVasicekproc.R
.draw.ouch: simulVasicekproc.R
.drawOneLineage: simulVasicekproc.R
.dzeta.matrix.t: precalcs.R
.dzeta.matrix.t1t2: precalcs.R
.dzeta.matrix.t1t2v2: precalcs.R
.dzeta.matrix.tv2: precalcs.R
.dzetaIJ.matrix: precalcs.R
.dzetaKappa.intercept: precalcs.R
.dzetaKappa.matrix: precalcs.R
.dzetaKappa.matrix.simp.t: precalcs.R
.dzetaKappa.matrix2: precalcs.R
.dzetaKappa.simp.intercept: precalcs.R
.EvaluatePoint: modelparams.R
.extract_from_signs: phylgls.R
.extract_GLS_results: phylgls.R
.extract_GLS_results_bm: phylgls.R
.extract_GLS_results_mvslouch: phylgls.R
.extract_GLS_results_ouch: phylgls.R
.extract.model.key.properties: evolmodelest.R
.generate.all.model.setups: evolmodelest.R
.generate.basic.model.setups: evolmodelest.R
.generate.ext.model.setups: evolmodelest.R
.generate.fund.model.setups: evolmodelest.R
.generate.list.of.model.setups: evolmodelest.R
.generate.model.description: evolmodelest.R
.generate.model.se: evolmodelest.R
.generate.univ.model.setups: evolmodelest.R
.generatenames: PhyloSDEestim.R
.getESScriteria: getESS.R
.getMVSLphylCovMatrix: getESS.R
.getMVSLproc_size: getESS.R
.get_fullGLSintercept: phylgls.R
.get_fullGLSintercept_mvslouch: phylgls.R
.givens.mult.M: matrixparametrizations.R
.givens.real: matrixparametrizations.R
.glsgc.estim: estimGLSGC.R
.InitialRegimeSetup: regimes.R
.internal_BrownianMotionModel: wrappers.R
.internal_mvslouchModel: wrappers.R
.internal_ouchModel: wrappers.R
.inv.rho.givens.real: matrixparametrizations.R
.is0_Merror: estimBM.R
.Kappa.matrix: precalcs.R
.Kappa.matrix.NA: precalcs.R
.M.mult.givens: matrixparametrizations.R
.make.states.mvsl: make.states.mvsl.R
.maxlik.estim: estimMAXLIK.R
.MinusPhylLogLikFunc: loglik.R
.MinusPhylLogLikFuncMVslouch: loglik.R
.MinusPhylLogLikFuncouch: loglik.R
.ModelSetupPCMBase: regimes.R
.mvslouch_to_ouch_model: modelparamstransform.R
.mvslouch.phyl.cov: sdecovariancephyl.R
.mvslouch.simulate: simulVasicekproc.R
.norm.max: modelparamssummary.R
.num.inv.transform.positive: matrixparametrizations.R
.num.transform.positive: matrixparametrizations.R
.oubm.sim.est: bootstrap.R
.ouch.phyl.cov: sdecovariancephyl.R
.ouch.simulate: simulVasicekproc.R
.ouch2slouch.mvsl: trees2slouch.mvsl.R
.ouou.sim.est: bootstrap.R
.OUtoBM_pcmbase_model_box: estimBM.R
.packageName:                              
.par_transform_withPCMBase: modelparamstransform.R
.par.inv.transform: modelparamstransform.R
.par.inv.transform.decomp.matrix: matrixparametrizations.R
.par.inv.transform.decomp.neg.real.matrix: matrixparametrizations.R
.par.inv.transform.decomp.pos.real.matrix: matrixparametrizations.R
.par.inv.transform.decomp.real.matrix: matrixparametrizations.R
.par.inv.transform.decomp.sym.matrix: matrixparametrizations.R
.par.inv.transform.invert.matrix: matrixparametrizations.R
.par.inv.transform.invert.qr.matrix matrixparametrizations.R
.par.inv.transform.invert.svd.matrix: matrixparametrizations.R
.par.inv.transform.lowertri.matrix: matrixparametrizations.R
.par.inv.transform.orth.matrix.cayley: matrixparametrizations.R
.par.inv.transform.orth.matrix.givens: matrixparametrizations.R
.par.inv.transform.symmetric.matrix: matrixparametrizations.R
.par.inv.transform.twobytwo.matrix: matrixparametrizations.R
.par.inv.transform.uppertri.matrix: matrixparametrizations.R
.par.transform: modelparamstransform.R
.par.transform.decomp.matrix: matrixparametrizations.R
.par.transform.decomp.neg.real.matrix: matrixparametrizations.R
.par.transform.decomp.pos.real.matrix: matrixparametrizations.R
.par.transform.decomp.real.matrix: matrixparametrizations.R
.par.transform.decomp.sym.matrix: matrixparametrizations.R
.par.transform.invert.matrix: matrixparametrizations.R
.par.transform.invert.qr.matrix: matrixparametrizations.R
.par.transform.invert.svd.matrix: matrixparametrizations.R
.par.transform.lowertri.matrix: matrixparametrizations.R
.par.transform.orth.matrix.cayley: matrixparametrizations.R
.par.transform.orth.matrix.givens: matrixparametrizations.R
.par.transform.symmetric.matrix: matrixparametrizations.R
.par.transform.twobytwo.matrix: matrixparametrizations.R
.par.transform.uppertri.matrix: matrixparametrizations.R
.Params.summary: modelparamssummary.R
.params.summary.bm: modelparamssummary.R
.params.summary.mvslouch: modelparamssummary.R
.params.summary.ouch: modelparamssummary.R
.pcmbaseDphylGaussian_RSS: phylgls.R
.pcmbaseDphylOU_GLS: phylgls.R
.PhyloSDEestim: PhyloSDEestim.R
.phyltree_remove_path_fields: phyltree_paths.R
.phyltree_remove_tips: phyltree_paths.R
.postproc_estres: evolmodelest.R
.return_best_model: evolmodelest.R
.row.rot: matrixparametrizations.R
.set_mean0_pcmbase_model_box: phygls.R
.set_1var_pcmbase_model_box: estimBM.R
.set.estimparams: PhyloSDEestim.R
.set.maxiter: PhyloSDEestim.R
.set.paramatrizationnames: PhyloSDEestim.R
.set_parscale: PhyloSDEestim.R
.set.tol: PhyloSDEestim.R
.simulVasicekProcPhylTree: simulVasicekprocphyl.R
.SummarizeFullPoint: modelparamssummary.R
.sym.par: matrixparametrizations.R
.sym.unpar: matrixparametrizations.R
.update_pcmbase_box_params: modelparamstransform.R
.update_pcmbase_box_params_bm: modelparamstransform.R
.update_pcmbase_box_params_mvslouch: modelparamstransform.R
.update_pcmbase_box_params_ouch: modelparamstransform.R
.internal_phyltree_paths_BL: phyltree_paths.R
.check_input_trait_data: PhyloSDEestim.R
.is_det0: matrixexps.R
.my_cov2cor: matrixexps.R
.my_stop: loglik.R
.my_message: loglik.R
.my_warning: loglik.R
.set_sing_bl_threshold: wrappers.R
.my_chol: evolmodelest.R
.set_pcmbase_model_box_X0: estimBM.R
.matrixcalc_is.diagonal.matrix: from_matrixcalc_1_0_3_matrixcalculations.R
.matrixcalc_is.positive.definite: from_matrixcalc_1_0_3_matrixcalculations.R
.matrixcalc_is.positive.semi.definite: from_matrixcalc_1_0_3_matrixcalculations.R
.matrixcalc_is.singular.matrix: from_matrixcalc_1_0_3_matrixcalculations.R
.matrixcalc_is.square.matrix: from_matrixcalc_1_0_3_matrixcalculations.R
.matrixcalc_is.symmetric.matrix: from_matrixcalc_1_0_3_matrixcalculations.R
.calc.phyl.mean: sdemeanphyl.R
.bm.phyl.mean: sdemeanphyl.R
.ouch.phyl.mean: sdemeanphyl.R
.mvslouch.phyl.mean: sdemeanphyl.R
.get_phyl_mean: OUphylregression.R
.OUphylreg_internal: OUphylregression.R
.id_evolmodel: OUphylregression.R
.f_warningsA: wrappers.R
}
\keyword{internal}
