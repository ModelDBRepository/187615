function MT = getMT_kharche()
 %       List of genes
 %         List of mutations
 %           List of groups of variables
 %             Pair (variable list + range)
 %               List of variables
 MT = []; 
 %CACNA1C:
 MT = [MT, { { { { 'offm_CaL12', -25.9 }, ...                                                  %http://www.ncbi.nlm.nih.gov/pubmed/19265197
                 { 'offh_CaL12', -27.0 } }, ...
               { { 'offm_CaL12', -37.3 }, ...                                                  %http://www.ncbi.nlm.nih.gov/pubmed/19265197
                 { 'offh_CaL12', -30.0 } },...
               { { 'offm_CaL12', [-31.4, +7.0] }, ...                                          %http://www.ncbi.nlm.nih.gov/pubmed/21685391
                 { 'slom_CaL12', [0.85, 1.45] }, ...                                           
                 { 'offh_CaL12', [-28.5, +16.3] }, ...                                         
                 { 'sloh_CaL12', [0.72, 1.38] } }, ...
               { { 'offm_CaL12', [-38.5, +12.9] }, ...                                         %http://www.ncbi.nlm.nih.gov/pubmed/16157588
                 { 'slom_CaL12', [0.46, 1.56] } },...
               { { 'offm_CaL12', [-27.8, +8.7] }, ...                                          %http://www.ncbi.nlm.nih.gov/pubmed/18836301
                 { 'slom_CaL12', [0.89, 1.14] },...
                 { 'offh_CaL12', [-19.1, +4.7] } }, ...                                          
               { { 'offm_CaL12', [-11.2, +1.0] }, ...                                          %http://www.ncbi.nlm.nih.gov/pubmed/15299022
                 { 'offh_CaL12', -3.1 },...
                 { 'sloh_CaL12', 1.24 } } } }];
 %CACNB2:
 MT = [MT, { { { { {'offh_CaL12', 'offh_CaL13'}, -5.2 }, ...                                   %http://www.ncbi.nlm.nih.gov/pubmed/19358333
                 { {'sloh_CaL12', 'sloh_CaL13'}, 0.69 } }, ... 
               { { 'tauh_CaL', 1.7} }, ...                                                     %http://www.ncbi.nlm.nih.gov/pubmed/7723731
               { { {'offm_CaL12', 'offm_CaL13'}, [-4.9, 4.9] }, ...                            %http://www.ncbi.nlm.nih.gov/pubmed/19723630
                 { {'offh_CaL12', 'offh_CaL13'}, [-5.1, 5.1] }, ...
                 { 'taum_CaL', [0.6, 1.68] }, ...
                 { 'tauh_CaL', [0.6, 1.66] } },...
               { { 'tauh_CaL', 1.26} } } }];                                                   %http://www.ncbi.nlm.nih.gov/pubmed/20025708
 %CACNA1D:
 MT = [MT, { { { { 'offm_CaL13', -10.9 }, ...                                                  %http://www.ncbi.nlm.nih.gov/pubmed/21998309 and
                 { 'slom_CaL13', 0.73 }, ...                                                   %http://www.ncbi.nlm.nih.gov/pubmed/21998310
                 { 'offh_CaL13', [-3.0, 3.5] }, ...                                            %(42A)
                 { 'sloh_CaL13', 0.81 }, ...
                 { 'tauh_CaL', 1.25 } }, ...
               { { 'offm_CaL13', [-10.6, 3.4] }, ...                                           %http://www.ncbi.nlm.nih.gov/pubmed/21998309 and
                 { 'slom_CaL13', [0.8, 1.12] }, ...                                            %http://www.ncbi.nlm.nih.gov/pubmed/21998310
                 { 'offh_CaL13', [-5.3, 1.2] }, ...                                            %(43S)
                 { 'sloh_CaL13', 0.66 }, ...
                 { 'tauh_CaL', 0.72 } }, ...
               { { 'offm_CaL13', 6.6 }, ...                                                    %http://www.ncbi.nlm.nih.gov/pubmed/20951705 and
                 { 'slom_CaL13', [0.75, 1.19] }, ...                                           %http://www.ncbi.nlm.nih.gov/pubmed/21054386
                 { 'tauh_CaL', [0.5, 1.12] } },...                                             %(CaV1.3 KO)
               { { 'offm_CaL13', -9.8 }, ...                                                   %http://www.ncbi.nlm.nih.gov/pubmed/25620733
                 { 'slom_CaL13', 0.8 }, ...                                        
                 { 'offh_CaL13', -15.4 }, ...                                        
                 { 'sloh_CaL13', 1.05 } },...                                          
               { { 'offm_CaL13', [-24.2, +6.1] }, ...                                          %http://www.ncbi.nlm.nih.gov/pubmed/23913004
                 { 'slom_CaL13', [0.7, 1.24] }, ...                                        
                 { 'offh_CaL13', -14.5 }, ...                                        
                 { 'sloh_CaL13', [0.72, 1.28] }, ...                                        
                 { 'tauh_CaL', 3.52 } },...                                          
               { { 'offm_CaL13', -17.8 }, ...                                                  %http://www.ncbi.nlm.nih.gov/pubmed/22760075
                 { 'slom_CaL13', 0.81 }, ...                                        
                 { 'tauh_CaL', [0.77, 1.31] } } } }];
 %CACNA1I:
 MT = [MT, { { { { 'offm_CaT', 1.3 }, ...                                                      %http://www.ncbi.nlm.nih.gov/pubmed/15254077
                 { 'offh_CaT', 1.6 }, ...
                 { 'taum_CaT', [0.87, 1.45] }, ...
                 { 'tauh_CaT', 0.8 } },...
               { { 'offm_CaT', -4.3 }, ...                                                     %http://www.ncbi.nlm.nih.gov/pubmed/12080115
                 { 'slom_CaT', 1.14 }, ...
                 { 'offh_CaT', -4.4 }, ...
                 { 'sloh_CaT', [0.89, 1.04] }, ...
                 { 'taum_CaT', 0.53 }, ...
                 { 'tauh_CaT', 0.46 } } } }];
 %CACNA1S:
 MT = [MT, { { { { 'taum_CaL', 0.67 } }, ...                                                   %http://www.ncbi.nlm.nih.gov/pubmed/20861472
               { { {'offm_CaL12', 'offm_CaL13'}, -30.02 }, ...                                 %http://www.ncbi.nlm.nih.gov/pubmed/19134469
                 { {'slom_CaL12', 'slom_CaL13'}, 0.62 }, ...
                 { 'taum_CaL', 0.49} } } }];
             
 %ATP2A:                                                                                                                                                                                                    
 MT = [MT, { { { { 'Pup_SERCA', 0.66 } },...                                                   %http://www.ncbi.nlm.nih.gov/pubmed/9891028
               { { 'Pup_SERCA', [0.22, 2.31] } },...                                           %http://www.ncbi.nlm.nih.gov/12975374
               { { 'Pup_SERCA', 0.2 } } } }];                                                  %http://www.ncbi.nlm.nih.gov/12670936

% %ATP2B:                                                                                                                                                                                                    
% MT.append([ [ [ 'decay_CaDynamics_E2', 1.97 ] ],                                             %http://www.ncbi.nlm.nih.gov/pubmed/22789621 (N/A)
%             [ [ 'decay_CaDynamics_E2', 1.5 ],                                                %http://www.ncbi.nlm.nih.gov/pubmed/21232211 (N/A)
%               [ 'minCai_CaDynamics_E2' , 1.4 ] ],
%             [ [ 'decay_CaDynamics_E2', 4.45 ] ] ])                                           %http://www.ncbi.nlm.nih.gov/pubmed/17234811 (N/A)

 %SCN1A:
 MT = [MT, { { { { 'offm_NaTTXS', -0.3 }, ...                                                  %http://www.ncbi.nlm.nih.gov/pubmed/18632931
                 { 'offh_NaTTXS', 5 }, ...
                 { 'slom_NaTTXS', 1.15 }, ...
                 { 'sloh_NaTTXS', 1.23 } }, ...
               { { 'offm_NaTTXS', 2.8 }, ...                                                   %http://www.ncbi.nlm.nih.gov/pubmed/17397047
                 { 'offh_NaTTXS', 9.6 }, ...
                 { 'slom_NaTTXS', 0.984 }, ...
                 { 'sloh_NaTTXS', 1.042 } }, ...
               { { 'offm_NaTTXS', -4.0 }, ...                                                  %http://www.ncbi.nlm.nih.gov/pubmed/21864321
                 { 'offh_NaTTXS', -5.8 }, ...
                 { 'slom_NaTTXS', 0.92 }, ...
                 { 'sloh_NaTTXS', 1.13 }, ...
                 { 'tauh_NaTTXS', 1.47 } }, ...
               { { 'offm_NaTTXS', -8.1 }, ...                                                  %http://www.ncbi.nlm.nih.gov/pubmed/21864321
                 { 'offh_NaTTXS', 2.2 }, ...
                 { 'slom_NaTTXS', 0.97 }, ...
                 { 'sloh_NaTTXS', 0.97 }, ...
                 { 'tauh_NaTTXS', 1.59 } }, ...
               { { 'offm_NaTTXS', 6.0 }, ...                                                   %http://www.ncbi.nlm.nih.gov/pubmed/23398611
                 { 'slom_NaTTXS', 1.16 }, ...
                 { 'tauh_NaTTXS', 1.29 } }, ...
               { { 'offm_NaTTXS', 10.0 }, ...                                                  %http://www.ncbi.nlm.nih.gov/pubmed/16326807
                 { 'offh_NaTTXS', -0.6 }, ...
                 { 'slom_NaTTXS', 1.15 }, ...
                 { 'sloh_NaTTXS', 1.14 } } } }];
%%SCN9A:
% MT = [MT, { { { { {'offh_NaTTXS'}, 0 } }, ...                                                 %http://www.ncbi.nlm.nih.gov/pubmed/22136189 (N/A)
%               { { {'offh_NaTTXS'}, 0 }, ...                                                   %http://www.ncbi.nlm.nih.gov/pubmed/18945915 (N/A)
%                 { 'sloh_NaTTXS', 1 }, ...
%                 { 'offm_NaTTXS', 0 }, ...
%                 { 'offh_NaTTXS', 0 }, ...
%                 { 'sloh_NaTTXS', 1 } }, ...
%               { { 'offm_NaTTXS', 0 }, ...                                                     %http://www.ncbi.nlm.nih.gov/pubmed/16392115 (N/A)
%                 { 'offh_NaTTXS', 0 } }, ...
%               { { 'offm_NaTTXS', 0 }, ...                                                     %http://www.ncbi.nlm.nih.gov/pubmed/15958509 (N/A)
%                 { 'offh_NaTTXS', 0 } } } }];
% %KCNS3:
% MT = [MT, { { { { {'taum_Kr'}, 1 }, ...                                                       %http://www.ncbi.nlm.nih.gov/pubmed/10484328 (N/A)
%                 { {'tauh1_Kr', 'tauh2_Kr'}, 1 }, ...
%                 { 'sloh_Kr', 1 } } } }];
% %KCNN3:
% MT = [MT, { { { { 'offc_SK_E2', 0.86 }, ...                                                   %http://www.ncbi.nlm.nih.gov/pubmed/14978258 (N/A)
%                { 'sloc_SK_E2', 1.24 } } } }];

 %HCN1:
 MT = [MT, { { { { 'offh_f', -26.5 }, ...                                                       %http://www.ncbi.nlm.nih.gov/pubmed/17185333
                 { 'sloh_f', 0.64 } }, ...
               { { 'offh_f', [-25.9, 17.7] }, ...                                               %http://www.ncbi.nlm.nih.gov/pubmed/12668666
                 { 'sloh_f', 0.6 } }, ...
               { { 'offh_f', 3.9 }, ...                                                         %http://www.ncbi.nlm.nih.gov/pubmed/26578877                                                                                                     
                 { 'tauh_f', 0.88 } } } }];

 %CACNB2 reprise:
 MT = [MT, { { { { {'offm_CaL12', 'offm_CaL13'}, 3 }, ...                            %http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N1 vs N4)
                 { {'offh_CaL12', 'offh_CaL13'}, 3.48 }, ...
                 { 'taum_CaL', 1.01 }, ...
                 { 'tauh_CaL', 0.89 } },...
               { { {'offm_CaL12', 'offm_CaL13'}, -1.11 }, ...                        %http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N3 vs N4)   
                 { {'offh_CaL12', 'offh_CaL13'}, 5.14 }, ...
                 { 'taum_CaL', 0.6 }, ...
                 { 'tauh_CaL', 1.48 } },...
               { { {'offm_CaL12', 'offm_CaL13'}, -1.88 }, ...                        %http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N5 vs N4)
                 { {'offh_CaL12', 'offh_CaL13'}, 2.69 }, ...
                 { 'taum_CaL', 0.6 }, ...
                 { 'tauh_CaL', 1.35 } } } }];
% %KCNB1:
% MT = [MT, { { { { 'offm_K_Pst', 5 }, ...                                                      %http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203K) (N/A)
%                 { 'offh_K_Pst', 3 }, ...
%                 { 'slom_K_Pst', 1.11 }, ...
%                 { 'sloh_K_Pst', 0.86 }, ...
%                 { {'taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'}, 0.5 }, ...
%                 { {'tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'}, 0.53 } }, ...
%               { { 'offm_K_Pst', 1 }, ...                                                      %http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203D) (N/A)
%                 { 'offh_K_Pst', -6 }, ...
%                 { 'slom_K_Pst', 1.22 }, ...
%                 { 'sloh_K_Pst', 1.0 }, ...
%                 { {'taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'}, 0.89 }, ...
%                 { {'tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'}, 1.13 } }, ...
%               { { 'offm_K_Pst', 6 }, ...                                                      %http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347K) (N/A)
%                 { 'offh_K_Pst', -8 }, ...
%                 { 'slom_K_Pst', 1.33 }, ...
%                 { 'sloh_K_Pst', 1.0 }, ...
%                 { {'taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'}, 0.5 }, ...
%                 { {'tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'}, 0.87 } }, ...
%               { { 'offm_K_Pst', -28 }, ...                                                    %http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347D) (N/A)
%                 { 'offh_K_Pst', -27 }, ...
%                 { 'slom_K_Pst', 1.11 }, ...
%                 { 'sloh_K_Pst', 0.71 }, ...
%                 { {'taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'}, 1.13 }, ...
%                 { {'tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'}, 2.27 } }, ... 
%               { { 'offm_K_Pst', 14 }, ...                                                     %http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203W) (N/A)
%                 { 'offh_K_Pst', -21 }, ...
%                 { 'slom_K_Pst', 2.0 }, ...
%                 { 'sloh_K_Pst', 1.0 }, ...
%                 { {'taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'}, 0.39 }, ...
%                 { {'tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'}, 1.2 } }, ... 
%               { { 'offm_K_Pst', -13 }, ...                                                    %http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347W) (N/A)
%                 { 'offh_K_Pst', -13 }, ...
%                 { 'slom_K_Pst', 1.33 }, ...
%                 { 'sloh_K_Pst', 0.71 }, ...
%                 { {'taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'}, 0.95 }, ...
%                 { {'tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'}, 5.13 } } } }];