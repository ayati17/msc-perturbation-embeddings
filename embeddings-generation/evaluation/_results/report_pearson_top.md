# Results report (pearson_top)

specs: 29 | models: 7 | rows: 48720


## Headline - absolute score and improvement over the mean baseline
                               spec   model  absolute       sd  n  ci95_abs   vs_mean  ci95_rel  beats_mean
        text_semantic/f2llm/curated      gp  0.441824 0.050516 30  0.018077  0.042341  0.008825        True
        text_semantic/qwen3/curated      gp  0.438989 0.055072 30  0.019707  0.039506  0.010050        True
        text_semantic/bekko/curated      gp  0.437309 0.058666 30  0.020993  0.037826  0.008821        True
         text_struct/bekko/combined      gp  0.432898 0.051350 30  0.018375  0.033415  0.008519        True
            text_struct/f2llm/iupac      gp  0.430621 0.063687 30  0.022790  0.031138  0.008697        True
                   string/chemberta      gp  0.430232 0.049256 30  0.017626  0.030749  0.009801        True
text_struct/f2llm/functional_groups      gp  0.429907 0.057271 30  0.020494  0.030424  0.007302        True
            text_struct/qwen3/iupac      gp  0.429795 0.061010 30  0.021832  0.030312  0.007441        True
text_struct/bekko/functional_groups      gp  0.428561 0.059200 30  0.021185  0.029078  0.006908        True
        text_semantic/bekko/curated   ridge  0.428503 0.055104 30  0.019719  0.029020  0.008786        True
        text_semantic/f2llm/curated   ridge  0.426753 0.049525 30  0.017722  0.027270  0.008568        True
         text_struct/f2llm/combined      gp  0.426688 0.059075 30  0.021140  0.027205  0.007463        True
        text_semantic/qwen3/curated   ridge  0.425949 0.053157 30  0.019022  0.026466  0.008577        True
         text_struct/qwen3/combined      gp  0.425005 0.059110 30  0.021152  0.025522  0.007814        True
                   string/selformer      gp  0.421761 0.058869 30  0.021066  0.022278  0.006460        True
                    ecfp_raw/counts      gp  0.421308 0.047516 30  0.017003  0.021826  0.011667        True
text_struct/qwen3/functional_groups      gp  0.417215 0.055390 30  0.019821  0.017732  0.005787        True
        text_semantic/bekko/curated     knn  0.416147 0.060687 30  0.021716  0.016664  0.013979        True
            text_struct/bekko/iupac      gp  0.415951 0.059328 30  0.021230  0.016468  0.007652        True
                             sand/z   ridge  0.414820 0.060320 30  0.021585  0.015338  0.009918        True
text_struct/bekko/functional_groups   ridge  0.414577 0.053507 30  0.019147  0.015094  0.006426        True
  text_struct/bekko/murcko_scaffold      gp  0.414305 0.056450 30  0.020200  0.014822  0.005567        True
                             sand/z      gp  0.412523 0.068319 30  0.024448  0.013040  0.014230       False
text_struct/f2llm/functional_groups   ridge  0.410862 0.057284 30  0.020499  0.011379  0.008017        True
  text_struct/qwen3/ecfp_positional      gp  0.410162 0.056343 30  0.020162  0.010679  0.006224        True
  text_struct/bekko/ecfp_positional      gp  0.410055 0.055236 30  0.019766  0.010572  0.006480        True
  text_struct/f2llm/ecfp_positional      gp  0.409718 0.059573 30  0.021318  0.010236  0.006255        True
  text_struct/f2llm/murcko_scaffold      gp  0.407690 0.056283 30  0.020141  0.008207  0.004074        True
  text_struct/qwen3/murcko_scaffold      gp  0.405775 0.055464 30  0.019848  0.006292  0.003944        True
        text_semantic/qwen3/curated     knn  0.405758 0.057278 30  0.020497  0.006275  0.012222       False
      text_struct/bekko/ecfp_binary      gp  0.405241 0.049889 30  0.017852  0.005758  0.006374       False
            text_struct/bekko/iupac   ridge  0.405239 0.053058 30  0.018986  0.005756  0.007193       False
         text_struct/bekko/combined   ridge  0.404703 0.050060 30  0.017914  0.005221  0.009577       False
        text_semantic/bekko/curated     mlp  0.403579 0.055876 30  0.019995  0.004096  0.014649       False
      text_struct/f2llm/ecfp_binary      gp  0.401758 0.052659 30  0.018844  0.002275  0.004937       False
            text_struct/qwen3/iupac   ridge  0.401670 0.056744 30  0.020305  0.002187  0.007883       False
                   string/selformer   ridge  0.400923 0.057448 30  0.020558  0.001441  0.008324       False
text_struct/bekko/rdkit_descriptors      gp  0.400805 0.057457 30  0.020561  0.001322  0.005390       False
                             sand/z     mlp  0.400381 0.066813 30  0.023909  0.000898  0.013891       False
text_struct/f2llm/rdkit_descriptors      gp  0.399784 0.056256 30  0.020131  0.000302  0.003319       False
text_struct/qwen3/functional_groups   ridge  0.399726 0.055864 30  0.019991  0.000243  0.008035       False
            text_struct/f2llm/iupac   ridge  0.398971 0.060316 30  0.021584 -0.000512  0.008349       False
         text_struct/f2llm/combined   ridge  0.398858 0.056237 30  0.020124 -0.000625  0.009174       False
  text_struct/bekko/ecfp_positional   ridge  0.398459 0.053086 30  0.018997 -0.001024  0.004995       False
      text_struct/qwen3/ecfp_binary      gp  0.398447 0.056438 30  0.020196 -0.001036  0.006127       False
text_struct/f2llm/rdkit_descriptors   ridge  0.397776 0.052579 30  0.018815 -0.001706  0.006056       False
        text_semantic/f2llm/curated     knn  0.396496 0.062266 30  0.022281 -0.002987  0.012763       False
                    ecfp_raw/counts     mlp  0.396106 0.049452 30  0.017696 -0.003377  0.013667       False
  text_struct/f2llm/ecfp_positional   ridge  0.396054 0.057961 30  0.020741 -0.003429  0.006375       False
         text_struct/qwen3/combined   ridge  0.396051 0.055948 30  0.020021 -0.003432  0.008956       False
text_struct/qwen3/rdkit_descriptors   ridge  0.395469 0.054349 30  0.019449 -0.004014  0.005867       False
                     string/chemgpt      gp  0.395348 0.049715 30  0.017790 -0.004135  0.005334       False
text_struct/qwen3/rdkit_descriptors      gp  0.394638 0.056477 30  0.020210 -0.004845  0.004364       False
text_struct/bekko/rdkit_descriptors   ridge  0.394607 0.054928 30  0.019656 -0.004876  0.004485       False
  text_struct/qwen3/murcko_scaffold   ridge  0.392618 0.056017 30  0.020046 -0.006865  0.006651       False
                    ecfp_raw/counts   ridge  0.392423 0.052711 30  0.018862 -0.007060  0.011350       False
                     string/chemgpt   ridge  0.391908 0.051299 30  0.018357 -0.007574  0.004377       False
      text_struct/f2llm/ecfp_binary   ridge  0.390783 0.052413 30  0.018756 -0.008700  0.005564       False
text_struct/f2llm/functional_groups     knn  0.390386 0.060087 30  0.021502 -0.009097  0.013960       False
  text_struct/bekko/murcko_scaffold   ridge  0.389992 0.055278 30  0.019781 -0.009491  0.005604       False
                   string/chemberta   ridge  0.389473 0.049732 30  0.017796 -0.010010  0.011131       False
  text_struct/f2llm/murcko_scaffold   ridge  0.389274 0.053183 30  0.019031 -0.010209  0.006242       False
         text_struct/bekko/combined     knn  0.386263 0.046668 30  0.016700 -0.013220  0.010166       False
        text_semantic/qwen3/curated     mlp  0.385404 0.058608 30  0.020973 -0.014079  0.016235       False
      text_struct/bekko/ecfp_binary   ridge  0.384111 0.051471 30  0.018419 -0.015371  0.005261       False
text_struct/bekko/functional_groups     knn  0.382694 0.054677 30  0.019566 -0.016789  0.011747       False
        text_semantic/f2llm/curated     mlp  0.381947 0.057244 30  0.020485 -0.017536  0.015539       False
  text_struct/qwen3/ecfp_positional   ridge  0.381803 0.055926 30  0.020013 -0.017680  0.007726       False
                   string/selformer     knn  0.381389 0.051536 30  0.018442 -0.018094  0.015447       False
                    ecfp_raw/counts     knn  0.380117 0.048854 30  0.017482 -0.019365  0.015340       False
         text_struct/f2llm/combined     knn  0.376101 0.056934 30  0.020374 -0.023382  0.011441       False
      text_struct/qwen3/ecfp_binary   ridge  0.374987 0.048335 30  0.017297 -0.024495  0.008004       False
                             sand/z     knn  0.374670 0.060445 30  0.021630 -0.024812  0.014367       False
            text_struct/f2llm/iupac     knn  0.374443 0.061351 30  0.021954 -0.025039  0.012950       False
                   string/chemberta     knn  0.371945 0.054119 30  0.019366 -0.027538  0.012477       False
         text_struct/qwen3/combined     knn  0.367559 0.057023 30  0.020406 -0.031924  0.012861       False
         text_struct/bekko/combined     mlp  0.365502 0.047369 30  0.016951 -0.033981  0.012361       False
            text_struct/bekko/iupac     knn  0.362604 0.050421 30  0.018043 -0.036879  0.010960       False
            text_struct/qwen3/iupac     knn  0.362262 0.055522 30  0.019868 -0.037221  0.012552       False
text_struct/qwen3/functional_groups     knn  0.359578 0.058224 30  0.020835 -0.039905  0.015297       False
  text_struct/qwen3/ecfp_positional     knn  0.359200 0.053632 30  0.019192 -0.040283  0.011489       False
  text_struct/f2llm/ecfp_positional     knn  0.358685 0.047375 30  0.016953 -0.040798  0.010348       False
            text_struct/bekko/iupac     mlp  0.353909 0.051469 30  0.018418 -0.045573  0.013344       False
  text_struct/bekko/ecfp_positional     knn  0.351740 0.053259 30  0.019058 -0.047743  0.012618       False
         text_struct/f2llm/combined     mlp  0.347681 0.062497 30  0.022364 -0.051801  0.015183       False
         text_struct/qwen3/combined     mlp  0.345833 0.067364 30  0.024106 -0.053650  0.017471       False
                   string/selformer     mlp  0.345332 0.063719 30  0.022801 -0.054151  0.016421       False
                   string/chemberta     mlp  0.345247 0.050730 30  0.018154 -0.054236  0.016127       False
            text_struct/f2llm/iupac     mlp  0.344116 0.063757 30  0.022815 -0.055367  0.014346       False
text_struct/bekko/functional_groups     mlp  0.343527 0.058572 30  0.020960 -0.055955  0.011302       False
            text_struct/qwen3/iupac     mlp  0.341587 0.058427 30  0.020908 -0.057896  0.014588       False
  text_struct/bekko/murcko_scaffold     knn  0.326701 0.063286 30  0.022647 -0.072782  0.014082       False
  text_struct/qwen3/ecfp_positional     mlp  0.325864 0.054939 30  0.019660 -0.073619  0.014921       False
  text_struct/f2llm/murcko_scaffold     knn  0.324534 0.055054 30  0.019701 -0.074949  0.011680       False
  text_struct/bekko/ecfp_positional     mlp  0.323950 0.046772 30  0.016737 -0.075533  0.012409       False
text_struct/f2llm/functional_groups     mlp  0.322271 0.065163 30  0.023318 -0.077212  0.024265       False
text_struct/bekko/rdkit_descriptors     mlp  0.319373 0.048143 30  0.017228 -0.080110  0.013852       False
  text_struct/qwen3/murcko_scaffold     knn  0.319041 0.051926 30  0.018582 -0.080442  0.013967       False
text_struct/bekko/rdkit_descriptors     knn  0.318998 0.047312 30  0.016930 -0.080485  0.012332       False
      text_struct/bekko/ecfp_binary     knn  0.318170 0.050226 30  0.017973 -0.081312  0.017336       False
  text_struct/f2llm/ecfp_positional     mlp  0.316456 0.051308 30  0.018360 -0.083027  0.016963       False
      text_struct/f2llm/ecfp_binary     knn  0.316401 0.053842 30  0.019267 -0.083081  0.016217       False
text_struct/f2llm/rdkit_descriptors     knn  0.314169 0.040263 30  0.014408 -0.085314  0.014154       False
                     string/chemgpt     knn  0.312994 0.050043 30  0.017908 -0.086489  0.012446       False
text_struct/qwen3/rdkit_descriptors     knn  0.311858 0.053001 30  0.018966 -0.087624  0.016409       False
      text_struct/qwen3/ecfp_binary     knn  0.310741 0.055007 30  0.019684 -0.088742  0.014394       False
text_struct/qwen3/functional_groups     mlp  0.303202 0.067000 30  0.023975 -0.096281  0.021031       False
  text_struct/bekko/murcko_scaffold     mlp  0.295058 0.062103 30  0.022223 -0.104425  0.019539       False
text_struct/f2llm/rdkit_descriptors     mlp  0.293010 0.045612 30  0.016322 -0.106473  0.015452       False
  text_struct/f2llm/murcko_scaffold     mlp  0.286446 0.068606 30  0.024550 -0.113036  0.017061       False
  text_struct/qwen3/murcko_scaffold     mlp  0.278865 0.055970 30  0.020029 -0.120618  0.015474       False
      text_struct/f2llm/ecfp_binary     mlp  0.277930 0.049431 30  0.017689 -0.121553  0.018610       False
      text_struct/qwen3/ecfp_binary     mlp  0.275218 0.046211 30  0.016536 -0.124265  0.017207       False
      text_struct/bekko/ecfp_binary     mlp  0.264356 0.047216 30  0.016896 -0.135127  0.020124       False
text_struct/qwen3/rdkit_descriptors     mlp  0.260021 0.049769 30  0.017810 -0.139462  0.017774       False
                     string/chemgpt     mlp  0.160202 0.099910 30  0.035752 -0.239281  0.040815       False
         text_struct/qwen3/combined  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
            text_struct/qwen3/iupac  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
text_struct/qwen3/rdkit_descriptors  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
  text_struct/qwen3/murcko_scaffold  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
  text_struct/qwen3/ecfp_positional  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
        text_semantic/bekko/curated  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
        text_semantic/qwen3/curated  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
        text_semantic/f2llm/curated  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
                   string/chemberta  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
                   string/selformer  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
                             sand/z  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
                     string/chemgpt  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
         text_struct/bekko/combined  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
      text_struct/bekko/ecfp_binary  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
                    ecfp_raw/counts  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
text_struct/f2llm/rdkit_descriptors  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
text_struct/qwen3/functional_groups  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
      text_struct/qwen3/ecfp_binary  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
  text_struct/f2llm/ecfp_positional  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
      text_struct/f2llm/ecfp_binary  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
text_struct/bekko/functional_groups  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
  text_struct/bekko/ecfp_positional  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
text_struct/bekko/rdkit_descriptors  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
            text_struct/bekko/iupac  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
text_struct/f2llm/functional_groups  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
            text_struct/f2llm/iupac  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
  text_struct/f2llm/murcko_scaffold  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
         text_struct/f2llm/combined  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
  text_struct/bekko/murcko_scaffold  random  0.150504 0.046295 30  0.016567 -0.248979  0.019817       False
        text_semantic/qwen3/curated control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
  text_struct/bekko/ecfp_positional control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
      text_struct/bekko/ecfp_binary control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
        text_semantic/bekko/curated control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
        text_semantic/f2llm/curated control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
                   string/chemberta control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
                             sand/z control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
                    ecfp_raw/counts control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
         text_struct/bekko/combined control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
                     string/chemgpt control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
                   string/selformer control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
text_struct/bekko/functional_groups control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
         text_struct/f2llm/combined control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
            text_struct/f2llm/iupac control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
text_struct/f2llm/functional_groups control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
      text_struct/f2llm/ecfp_binary control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
  text_struct/f2llm/ecfp_positional control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
            text_struct/bekko/iupac control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
  text_struct/bekko/murcko_scaffold control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
text_struct/bekko/rdkit_descriptors control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
  text_struct/f2llm/murcko_scaffold control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
  text_struct/qwen3/ecfp_positional control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
      text_struct/qwen3/ecfp_binary control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
text_struct/f2llm/rdkit_descriptors control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
         text_struct/qwen3/combined control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
            text_struct/qwen3/iupac control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
text_struct/qwen3/functional_groups control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
  text_struct/qwen3/murcko_scaffold control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False
text_struct/qwen3/rdkit_descriptors control  0.000000 0.000000 30  0.000000 -0.399483  0.019196       False


## All metrics (mean over folds and cell lines)
metric                                       direction    mae    mse  overlap_top  pearson  pearson_top      r2  r2_gene_sd
spec                                model                                                                                  
ecfp_raw/counts                     control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.619  0.022  0.004        0.429    0.331        0.421  -0.499       0.501
                                    knn          0.646  0.023  0.004        0.410    0.295        0.380  -0.909       0.629
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.610  0.022  0.004        0.430    0.312        0.396  -0.962       0.678
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.623  0.022  0.004        0.423    0.308        0.392  -0.634       0.444
sand/z                              control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.619  0.021  0.004        0.427    0.327        0.413  -0.463       0.570
                                    knn          0.649  0.023  0.004        0.413    0.292        0.375  -0.875       0.683
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.610  0.022  0.004        0.424    0.315        0.400  -0.505       0.598
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.624  0.022  0.004        0.426    0.326        0.415  -0.454       0.362
string/chemberta                    control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.021  0.004        0.433    0.339        0.430  -0.407       0.457
                                    knn          0.650  0.023  0.004        0.419    0.291        0.372  -0.847       0.673
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.595  0.024  0.005        0.420    0.272        0.345  -1.233       0.788
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.623  0.022  0.004        0.422    0.307        0.389  -0.552       0.392
string/chemgpt                      control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.619  0.022  0.004        0.411    0.301        0.395  -0.366      -0.152
                                    knn          0.642  0.024  0.005        0.397    0.231        0.313  -1.133       0.587
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.508  0.185  0.070        0.143    0.080        0.160 -51.585     -14.596
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.618  0.022  0.005        0.409    0.298        0.392  -0.401      -0.097
string/selformer                    control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.621  0.022  0.004        0.427    0.330        0.422  -0.436       0.214
                                    knn          0.649  0.023  0.004        0.414    0.299        0.381  -0.674       0.607
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.584  0.025  0.005        0.412    0.267        0.345  -1.635       0.808
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.623  0.022  0.004        0.421    0.317        0.401  -0.589       0.342
text_semantic/bekko/curated         control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.625  0.021  0.004        0.441    0.348        0.437  -0.978       0.641
                                    knn          0.656  0.022  0.004        0.432    0.332        0.416  -1.008       0.757
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.608  0.022  0.004        0.427    0.323        0.404  -0.933       0.691
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.628  0.021  0.004        0.432    0.345        0.429  -0.598       0.425
text_semantic/f2llm/curated         control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.623  0.021  0.004        0.438    0.355        0.442  -0.466       0.537
                                    knn          0.652  0.022  0.004        0.424    0.318        0.396  -0.902       0.718
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.596  0.027  0.005        0.412    0.304        0.382  -1.802       0.297
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.627  0.021  0.004        0.432    0.343        0.427  -0.487       0.430
text_semantic/qwen3/curated         control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.021  0.004        0.439    0.352        0.439  -0.745       0.584
                                    knn          0.653  0.022  0.004        0.427    0.324        0.406  -0.990       0.744
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.596  0.028  0.006        0.413    0.303        0.385  -2.127       0.446
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.626  0.021  0.004        0.429    0.340        0.426  -0.529       0.381
text_struct/bekko/combined          control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.625  0.021  0.004        0.438    0.343        0.433  -0.491       0.490
                                    knn          0.650  0.023  0.004        0.427    0.301        0.386  -1.152       0.722
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.598  0.023  0.004        0.423    0.288        0.366  -1.098       0.742
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.625  0.022  0.004        0.428    0.321        0.405  -0.549       0.385
text_struct/bekko/ecfp_binary       control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.624  0.022  0.004        0.411    0.309        0.405  -0.309       0.116
                                    knn          0.644  0.024  0.005        0.394    0.239        0.318  -0.980       0.601
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.583  0.027  0.007        0.394    0.195        0.264  -2.350       0.799
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.619  0.022  0.005        0.406    0.292        0.384  -0.416       0.003
text_struct/bekko/ecfp_positional   control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.621  0.022  0.004        0.422    0.318        0.410  -0.405       0.188
                                    knn          0.647  0.023  0.005        0.405    0.272        0.352  -0.949       0.656
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.595  0.024  0.005        0.409    0.250        0.324  -1.273       0.707
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.620  0.022  0.004        0.418    0.310        0.398  -0.410       0.130
text_struct/bekko/functional_groups control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.021  0.004        0.438    0.342        0.429  -0.516       0.385
                                    knn          0.652  0.022  0.004        0.424    0.302        0.383  -0.911       0.710
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.588  0.025  0.005        0.418    0.272        0.344  -2.031       0.830
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.625  0.021  0.004        0.431    0.330        0.415  -0.613       0.381
text_struct/bekko/iupac             control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.620  0.022  0.004        0.431    0.327        0.416  -0.593       0.357
                                    knn          0.648  0.023  0.005        0.418    0.279        0.363  -1.447       0.692
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.596  0.024  0.005        0.420    0.275        0.354  -1.682       0.770
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.623  0.022  0.004        0.425    0.317        0.405  -0.639       0.287
text_struct/bekko/murcko_scaffold   control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.022  0.004        0.420    0.320        0.414  -0.376       0.033
                                    knn          0.643  0.025  0.005        0.401    0.246        0.327  -1.197       0.628
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.585  0.027  0.006        0.402    0.223        0.295  -2.617       0.805
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.619  0.022  0.005        0.412    0.296        0.390  -0.644       0.066
text_struct/bekko/rdkit_descriptors control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.618  0.022  0.004        0.413    0.308        0.401  -0.339       0.117
                                    knn          0.643  0.024  0.005        0.397    0.240        0.319  -0.822       0.610
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.586  0.025  0.005        0.398    0.233        0.319  -1.103       0.675
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.620  0.022  0.004        0.412    0.304        0.395  -0.336       0.034
text_struct/f2llm/combined          control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.623  0.022  0.004        0.436    0.337        0.427  -0.505       0.416
                                    knn          0.649  0.023  0.004        0.426    0.296        0.376  -1.268       0.701
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.585  0.024  0.004        0.416    0.270        0.348  -1.329       0.821
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.624  0.022  0.004        0.425    0.314        0.399  -0.603       0.386
text_struct/f2llm/ecfp_binary       control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.623  0.022  0.004        0.415    0.308        0.402  -0.347       0.038
                                    knn          0.643  0.024  0.005        0.393    0.237        0.316  -0.868       0.561
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.573  0.026  0.006        0.393    0.212        0.278  -1.783       0.783
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.618  0.022  0.004        0.408    0.298        0.391  -0.366      -0.005
text_struct/f2llm/ecfp_positional   control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.022  0.004        0.425    0.317        0.410  -0.464       0.258
                                    knn          0.648  0.023  0.004        0.409    0.275        0.359  -0.913       0.662
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.579  0.025  0.005        0.405    0.245        0.316  -1.375       0.768
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.621  0.022  0.004        0.417    0.308        0.396  -0.443       0.186
text_struct/f2llm/functional_groups control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.021  0.004        0.436    0.342        0.430  -0.459       0.419
                                    knn          0.652  0.022  0.004        0.428    0.312        0.390  -0.894       0.729
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.566  0.037  0.008        0.380    0.243        0.322  -5.642       0.357
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.625  0.021  0.004        0.430    0.326        0.411  -0.576       0.407
text_struct/f2llm/iupac             control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.623  0.022  0.004        0.435    0.339        0.431  -0.490       0.379
                                    knn          0.647  0.023  0.005        0.419    0.289        0.374  -1.328       0.689
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.585  0.024  0.005        0.415    0.266        0.344  -1.540       0.818
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.623  0.022  0.004        0.420    0.311        0.399  -0.567       0.330
text_struct/f2llm/murcko_scaffold   control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.022  0.004        0.414    0.313        0.408  -0.347      -0.070
                                    knn          0.639  0.025  0.006        0.398    0.240        0.325  -1.765       0.628
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.574  0.028  0.006        0.396    0.213        0.286  -2.700       0.831
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.619  0.022  0.005        0.409    0.295        0.389  -0.525       0.092
text_struct/f2llm/rdkit_descriptors control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.620  0.022  0.004        0.412    0.306        0.400  -0.377      -0.055
                                    knn          0.643  0.024  0.005        0.395    0.235        0.314  -0.763       0.605
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.577  0.026  0.006        0.391    0.217        0.293  -1.484       0.762
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.621  0.022  0.004        0.409    0.307        0.398  -0.355       0.030
text_struct/qwen3/combined          control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.621  0.022  0.004        0.428    0.334        0.425  -0.571       0.301
                                    knn          0.647  0.023  0.005        0.417    0.285        0.368  -1.253       0.667
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.584  0.028  0.006        0.402    0.268        0.346  -2.040       0.471
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.622  0.022  0.004        0.421    0.311        0.396  -0.611       0.316
text_struct/qwen3/ecfp_binary       control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.022  0.004        0.411    0.303        0.398  -0.337       0.019
                                    knn          0.642  0.024  0.005        0.392    0.231        0.311  -0.796       0.589
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.575  0.026  0.006        0.393    0.211        0.275  -1.661       0.816
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.620  0.022  0.004        0.408    0.286        0.375  -0.448       0.107
text_struct/qwen3/ecfp_positional   control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.621  0.022  0.004        0.419    0.316        0.410  -0.394       0.175
                                    knn          0.646  0.023  0.005        0.409    0.272        0.359  -0.907       0.639
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.577  0.025  0.005        0.404    0.248        0.326  -1.263       0.765
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.622  0.022  0.004        0.414    0.299        0.382  -0.501       0.232
text_struct/qwen3/functional_groups control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.022  0.004        0.422    0.323        0.417  -0.418       0.116
                                    knn          0.645  0.024  0.005        0.410    0.274        0.360  -1.183       0.630
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.572  0.027  0.006        0.403    0.231        0.303  -2.117       0.858
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.620  0.022  0.004        0.419    0.309        0.400  -0.467       0.109
text_struct/qwen3/iupac             control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.623  0.022  0.004        0.433    0.338        0.430  -0.513       0.371
                                    knn          0.647  0.023  0.004        0.416    0.281        0.362  -1.151       0.688
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.587  0.024  0.005        0.414    0.264        0.342  -1.493       0.789
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.623  0.022  0.004        0.422    0.312        0.402  -0.586       0.321
text_struct/qwen3/murcko_scaffold   control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.622  0.022  0.004        0.414    0.311        0.406  -0.338      -0.116
                                    knn          0.641  0.025  0.005        0.392    0.234        0.319  -1.317       0.617
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.580  0.027  0.007        0.393    0.207        0.279  -2.976       0.808
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.619  0.022  0.005        0.411    0.298        0.393  -0.474       0.006
text_struct/qwen3/rdkit_descriptors control      0.044  0.022  0.005        0.133    0.000        0.000  -0.012      -0.284
                                    gp           0.620  0.022  0.004        0.411    0.302        0.395  -0.400      -0.058
                                    knn          0.643  0.024  0.005        0.399    0.233        0.312  -0.938       0.624
                                    mean         0.620  0.022  0.004        0.411    0.305        0.399  -0.336      -0.284
                                    mlp          0.574  0.027  0.006        0.389    0.192        0.260  -1.874       0.819
                                    random       0.688  0.030  0.008        0.358    0.111        0.151  -3.824       0.866
                                    ridge        0.620  0.022  0.004        0.410    0.305        0.395  -0.373       0.010


## By cell line 
cell_line                                     A549   K562   MCF7
spec                                model                       
ecfp_raw/counts                     control  0.000  0.000  0.000
                                    gp       0.417  0.407  0.440
                                    knn      0.372  0.357  0.412
                                    mean     0.394  0.381  0.423
                                    mlp      0.394  0.381  0.414
                                    random   0.147  0.130  0.175
                                    ridge    0.389  0.372  0.416
sand/z                              control  0.000  0.000  0.000
                                    gp       0.421  0.384  0.433
                                    knn      0.372  0.341  0.411
                                    mean     0.394  0.381  0.423
                                    mlp      0.405  0.363  0.433
                                    random   0.147  0.130  0.175
                                    ridge    0.421  0.385  0.438
string/chemberta                    control  0.000  0.000  0.000
                                    gp       0.421  0.418  0.451
                                    knn      0.357  0.352  0.407
                                    mean     0.394  0.381  0.423
                                    mlp      0.332  0.332  0.372
                                    random   0.147  0.130  0.175
                                    ridge    0.383  0.371  0.414
string/chemgpt                      control  0.000  0.000  0.000
                                    gp       0.387  0.382  0.417
                                    knn      0.296  0.292  0.351
                                    mean     0.394  0.381  0.423
                                    mlp      0.135  0.125  0.221
                                    random   0.147  0.130  0.175
                                    ridge    0.380  0.378  0.418
string/selformer                    control  0.000  0.000  0.000
                                    gp       0.413  0.396  0.456
                                    knn      0.376  0.338  0.430
                                    mean     0.394  0.381  0.423
                                    mlp      0.335  0.322  0.380
                                    random   0.147  0.130  0.175
                                    ridge    0.389  0.387  0.426
text_semantic/bekko/curated         control  0.000  0.000  0.000
                                    gp       0.439  0.420  0.453
                                    knn      0.417  0.393  0.438
                                    mean     0.394  0.381  0.423
                                    mlp      0.405  0.375  0.430
                                    random   0.147  0.130  0.175
                                    ridge    0.432  0.409  0.445
text_semantic/f2llm/curated         control  0.000  0.000  0.000
                                    gp       0.443  0.419  0.464
                                    knn      0.395  0.365  0.429
                                    mean     0.394  0.381  0.423
                                    mlp      0.364  0.355  0.426
                                    random   0.147  0.130  0.175
                                    ridge    0.426  0.410  0.444
text_semantic/qwen3/curated         control  0.000  0.000  0.000
                                    gp       0.446  0.410  0.461
                                    knn      0.407  0.377  0.434
                                    mean     0.394  0.381  0.423
                                    mlp      0.374  0.361  0.421
                                    random   0.147  0.130  0.175
                                    ridge    0.430  0.403  0.445
text_struct/bekko/combined          control  0.000  0.000  0.000
                                    gp       0.431  0.417  0.451
                                    knn      0.382  0.373  0.404
                                    mean     0.394  0.381  0.423
                                    mlp      0.369  0.335  0.392
                                    random   0.147  0.130  0.175
                                    ridge    0.406  0.388  0.420
text_struct/bekko/ecfp_binary       control  0.000  0.000  0.000
                                    gp       0.403  0.392  0.421
                                    knn      0.291  0.306  0.357
                                    mean     0.394  0.381  0.423
                                    mlp      0.265  0.258  0.270
                                    random   0.147  0.130  0.175
                                    ridge    0.375  0.369  0.408
text_struct/bekko/ecfp_positional   control  0.000  0.000  0.000
                                    gp       0.407  0.391  0.432
                                    knn      0.333  0.342  0.380
                                    mean     0.394  0.381  0.423
                                    mlp      0.318  0.320  0.334
                                    random   0.147  0.130  0.175
                                    ridge    0.385  0.388  0.422
text_struct/bekko/functional_groups control  0.000  0.000  0.000
                                    gp       0.424  0.410  0.453
                                    knn      0.371  0.365  0.411
                                    mean     0.394  0.381  0.423
                                    mlp      0.342  0.323  0.366
                                    random   0.147  0.130  0.175
                                    ridge    0.411  0.395  0.438
text_struct/bekko/iupac             control  0.000  0.000  0.000
                                    gp       0.405  0.396  0.447
                                    knn      0.359  0.337  0.392
                                    mean     0.394  0.381  0.423
                                    mlp      0.352  0.328  0.382
                                    random   0.147  0.130  0.175
                                    ridge    0.402  0.386  0.428
text_struct/bekko/murcko_scaffold   control  0.000  0.000  0.000
                                    gp       0.410  0.390  0.443
                                    knn      0.324  0.286  0.371
                                    mean     0.394  0.381  0.423
                                    mlp      0.285  0.275  0.325
                                    random   0.147  0.130  0.175
                                    ridge    0.379  0.375  0.416
text_struct/bekko/rdkit_descriptors control  0.000  0.000  0.000
                                    gp       0.395  0.384  0.424
                                    knn      0.304  0.295  0.358
                                    mean     0.394  0.381  0.423
                                    mlp      0.299  0.311  0.348
                                    random   0.147  0.130  0.175
                                    ridge    0.387  0.376  0.422
text_struct/f2llm/combined          control  0.000  0.000  0.000
                                    gp       0.425  0.405  0.450
                                    knn      0.362  0.353  0.413
                                    mean     0.394  0.381  0.423
                                    mlp      0.345  0.313  0.386
                                    random   0.147  0.130  0.175
                                    ridge    0.401  0.377  0.419
text_struct/f2llm/ecfp_binary       control  0.000  0.000  0.000
                                    gp       0.397  0.386  0.422
                                    knn      0.300  0.302  0.347
                                    mean     0.394  0.381  0.423
                                    mlp      0.277  0.249  0.308
                                    random   0.147  0.130  0.175
                                    ridge    0.378  0.374  0.421
text_struct/f2llm/ecfp_positional   control  0.000  0.000  0.000
                                    gp       0.405  0.394  0.430
                                    knn      0.349  0.341  0.387
                                    mean     0.394  0.381  0.423
                                    mlp      0.294  0.307  0.348
                                    random   0.147  0.130  0.175
                                    ridge    0.392  0.385  0.412
text_struct/f2llm/functional_groups control  0.000  0.000  0.000
                                    gp       0.427  0.411  0.451
                                    knn      0.377  0.372  0.422
                                    mean     0.394  0.381  0.423
                                    mlp      0.303  0.309  0.355
                                    random   0.147  0.130  0.175
                                    ridge    0.411  0.390  0.431
text_struct/f2llm/iupac             control  0.000  0.000  0.000
                                    gp       0.427  0.403  0.461
                                    knn      0.367  0.336  0.420
                                    mean     0.394  0.381  0.423
                                    mlp      0.339  0.310  0.383
                                    random   0.147  0.130  0.175
                                    ridge    0.397  0.374  0.426
text_struct/f2llm/murcko_scaffold   control  0.000  0.000  0.000
                                    gp       0.406  0.388  0.429
                                    knn      0.313  0.299  0.362
                                    mean     0.394  0.381  0.423
                                    mlp      0.265  0.267  0.327
                                    random   0.147  0.130  0.175
                                    ridge    0.379  0.376  0.413
text_struct/f2llm/rdkit_descriptors control  0.000  0.000  0.000
                                    gp       0.396  0.384  0.419
                                    knn      0.292  0.309  0.342
                                    mean     0.394  0.381  0.423
                                    mlp      0.271  0.282  0.326
                                    random   0.147  0.130  0.175
                                    ridge    0.391  0.384  0.418
text_struct/qwen3/combined          control  0.000  0.000  0.000
                                    gp       0.421  0.401  0.453
                                    knn      0.361  0.342  0.400
                                    mean     0.394  0.381  0.423
                                    mlp      0.340  0.309  0.389
                                    random   0.147  0.130  0.175
                                    ridge    0.397  0.372  0.420
text_struct/qwen3/ecfp_binary       control  0.000  0.000  0.000
                                    gp       0.400  0.378  0.417
                                    knn      0.296  0.288  0.348
                                    mean     0.394  0.381  0.423
                                    mlp      0.278  0.262  0.285
                                    random   0.147  0.130  0.175
                                    ridge    0.360  0.357  0.409
text_struct/qwen3/ecfp_positional   control  0.000  0.000  0.000
                                    gp       0.408  0.393  0.430
                                    knn      0.343  0.348  0.387
                                    mean     0.394  0.381  0.423
                                    mlp      0.314  0.314  0.350
                                    random   0.147  0.130  0.175
                                    ridge    0.382  0.365  0.398
text_struct/qwen3/functional_groups control  0.000  0.000  0.000
                                    gp       0.410  0.395  0.446
                                    knn      0.341  0.346  0.391
                                    mean     0.394  0.381  0.423
                                    mlp      0.285  0.280  0.345
                                    random   0.147  0.130  0.175
                                    ridge    0.394  0.380  0.424
text_struct/qwen3/iupac             control  0.000  0.000  0.000
                                    gp       0.428  0.403  0.458
                                    knn      0.359  0.332  0.396
                                    mean     0.394  0.381  0.423
                                    mlp      0.341  0.313  0.371
                                    random   0.147  0.130  0.175
                                    ridge    0.403  0.377  0.426
text_struct/qwen3/murcko_scaffold   control  0.000  0.000  0.000
                                    gp       0.404  0.388  0.426
                                    knn      0.313  0.292  0.352
                                    mean     0.394  0.381  0.423
                                    mlp      0.256  0.266  0.315
                                    random   0.147  0.130  0.175
                                    ridge    0.376  0.381  0.421
text_struct/qwen3/rdkit_descriptors control  0.000  0.000  0.000
                                    gp       0.390  0.378  0.416
                                    knn      0.274  0.307  0.355
                                    mean     0.394  0.381  0.423
                                    mlp      0.241  0.252  0.287
                                    random   0.147  0.130  0.175
                                    ridge    0.385  0.383  0.418


## random vs scaffold split 
split                                        random  scaffold    gap
spec                                model                           
string/chemgpt                      mlp       0.182     0.138  0.044
sand/z                              mlp       0.419     0.381  0.038
text_struct/f2llm/ecfp_binary       mlp       0.296     0.260  0.036
text_struct/qwen3/murcko_scaffold   mlp       0.296     0.261  0.035
text_struct/f2llm/ecfp_binary       knn       0.334     0.299  0.035
text_struct/f2llm/murcko_scaffold   knn       0.342     0.307  0.034
text_struct/bekko/ecfp_binary       mlp       0.281     0.247  0.034
text_struct/bekko/murcko_scaffold   mlp       0.311     0.279  0.031
text_struct/bekko/ecfp_binary       knn       0.333     0.304  0.029
text_struct/bekko/murcko_scaffold   knn       0.341     0.313  0.028
sand/z                              gp        0.426     0.399  0.027
                                    ridge     0.428     0.401  0.027
string/chemberta                    knn       0.385     0.359  0.026
text_struct/qwen3/murcko_scaffold   knn       0.331     0.307  0.025
string/selformer                    mlp       0.358     0.333  0.025
ecfp_raw/counts                     mlp       0.408     0.384  0.024
text_struct/qwen3/iupac             knn       0.374     0.350  0.024
text_struct/qwen3/functional_groups knn       0.372     0.348  0.024
text_struct/f2llm/rdkit_descriptors knn       0.326     0.302  0.023
text_struct/f2llm/ecfp_positional   mlp       0.328     0.305  0.023
text_struct/bekko/functional_groups mlp       0.355     0.332  0.023
text_struct/qwen3/ecfp_binary       mlp       0.286     0.264  0.022
text_struct/f2llm/rdkit_descriptors mlp       0.304     0.282  0.022
text_struct/bekko/ecfp_positional   knn       0.362     0.341  0.021
text_struct/qwen3/functional_groups mlp       0.314     0.293  0.021
text_struct/qwen3/ecfp_positional   gp        0.421     0.400  0.021
text_struct/qwen3/rdkit_descriptors knn       0.322     0.302  0.020
text_struct/f2llm/ecfp_binary       gp        0.411     0.392  0.019
text_struct/bekko/ecfp_positional   gp        0.419     0.401  0.018
text_struct/f2llm/ecfp_positional   gp        0.419     0.401  0.018
text_struct/f2llm/murcko_scaffold   ridge     0.398     0.380  0.018
                                    mlp       0.295     0.278  0.018
string/chemberta                    gp        0.439     0.421  0.018
text_struct/qwen3/murcko_scaffold   ridge     0.401     0.384  0.018
ecfp_raw/counts                     gp        0.430     0.413  0.017
text_struct/bekko/murcko_scaffold   ridge     0.399     0.381  0.017
text_struct/bekko/ecfp_binary       gp        0.414     0.397  0.017
text_struct/bekko/murcko_scaffold   gp        0.423     0.406  0.017
ecfp_raw/counts                     ridge     0.401     0.384  0.017
text_struct/qwen3/murcko_scaffold   gp        0.414     0.397  0.017
text_struct/f2llm/ecfp_binary       ridge     0.399     0.383  0.016
text_struct/bekko/rdkit_descriptors knn       0.327     0.311  0.016
text_struct/f2llm/iupac             gp        0.439     0.423  0.016
text_struct/bekko/functional_groups knn       0.391     0.375  0.016
text_struct/bekko/iupac             mlp       0.362     0.346  0.016
text_struct/f2llm/iupac             knn       0.382     0.366  0.016
sand/z                              knn       0.382     0.367  0.015
text_struct/qwen3/rdkit_descriptors gp        0.402     0.387  0.015
text_struct/f2llm/murcko_scaffold   gp        0.415     0.400  0.015
text_struct/f2llm/iupac             ridge     0.406     0.392  0.015
text_struct/qwen3/functional_groups gp        0.424     0.410  0.015
string/chemberta                    ridge     0.396     0.383  0.014
ecfp_raw/counts                     knn       0.387     0.373  0.013
text_struct/f2llm/combined          gp        0.433     0.420  0.013
text_struct/qwen3/combined          ridge     0.403     0.390  0.013
text_struct/f2llm/functional_groups gp        0.436     0.423  0.013
text_struct/qwen3/iupac             gp        0.436     0.423  0.013
text_struct/f2llm/ecfp_positional   ridge     0.402     0.390  0.013
text_struct/qwen3/ecfp_binary       gp        0.405     0.392  0.012
text_struct/qwen3/combined          gp        0.431     0.419  0.012
text_struct/bekko/iupac             gp        0.422     0.410  0.012
text_struct/qwen3/ecfp_binary       knn       0.317     0.305  0.012
text_struct/bekko/ecfp_binary       ridge     0.390     0.378  0.012
text_struct/qwen3/ecfp_positional   knn       0.365     0.353  0.012
text_struct/f2llm/functional_groups knn       0.396     0.384  0.012
text_struct/f2llm/rdkit_descriptors gp        0.406     0.394  0.012
text_struct/bekko/functional_groups gp        0.434     0.423  0.011
text_struct/qwen3/combined          mlp       0.351     0.340  0.011
text_struct/f2llm/functional_groups ridge     0.416     0.405  0.011
text_struct/f2llm/combined          knn       0.382     0.371  0.011
text_semantic/bekko/curated         gp        0.443     0.432  0.011
text_struct/bekko/ecfp_positional   ridge     0.404     0.393  0.010
text_struct/qwen3/functional_groups ridge     0.405     0.395  0.010
text_struct/bekko/iupac             ridge     0.410     0.400  0.010
text_struct/bekko/rdkit_descriptors gp        0.406     0.396  0.010
text_semantic/bekko/curated         mlp       0.409     0.399  0.010
string/selformer                    knn       0.386     0.377  0.009
text_struct/qwen3/ecfp_positional   ridge     0.386     0.377  0.009
text_semantic/qwen3/curated         mlp       0.390     0.381  0.009
string/chemberta                    mlp       0.350     0.341  0.009
text_struct/f2llm/ecfp_positional   knn       0.363     0.354  0.009
text_semantic/f2llm/curated         gp        0.446     0.437  0.009
text_struct/bekko/combined          ridge     0.409     0.400  0.009
text_struct/bekko/functional_groups ridge     0.419     0.410  0.009
text_semantic/qwen3/curated         ridge     0.430     0.422  0.009
text_struct/f2llm/combined          ridge     0.403     0.395  0.008
text_struct/bekko/rdkit_descriptors ridge     0.399     0.390  0.008
text_semantic/qwen3/curated         gp        0.443     0.435  0.008
text_struct/bekko/ecfp_positional   mlp       0.328     0.320  0.008
text_struct/bekko/iupac             knn       0.367     0.359  0.008
text_struct/qwen3/ecfp_positional   mlp       0.330     0.322  0.008
string/selformer                    ridge     0.405     0.397  0.008
text_semantic/qwen3/curated         knn       0.410     0.402  0.008
string/selformer                    gp        0.425     0.418  0.007
text_struct/qwen3/iupac             ridge     0.405     0.398  0.007
text_struct/bekko/combined          mlp       0.369     0.362  0.007
text_struct/qwen3/rdkit_descriptors ridge     0.399     0.392  0.006
text_struct/bekko/rdkit_descriptors mlp       0.322     0.317  0.006
text_semantic/f2llm/curated         ridge     0.429     0.424  0.005
text_struct/qwen3/ecfp_binary       ridge     0.378     0.372  0.005
text_struct/qwen3/rdkit_descriptors mlp       0.262     0.258  0.004
string/chemgpt                      ridge     0.394     0.390  0.004
text_struct/bekko/combined          gp        0.435     0.431  0.004
                                    knn       0.388     0.384  0.004
text_struct/qwen3/ecfp_positional   mean      0.401     0.398  0.003
text_struct/qwen3/functional_groups mean      0.401     0.398  0.003
text_struct/qwen3/murcko_scaffold   mean      0.401     0.398  0.003
string/chemgpt                      mean      0.401     0.398  0.003
text_struct/bekko/combined          mean      0.401     0.398  0.003
text_struct/bekko/ecfp_positional   mean      0.401     0.398  0.003
sand/z                              mean      0.401     0.398  0.003
ecfp_raw/counts                     mean      0.401     0.398  0.003
string/chemberta                    mean      0.401     0.398  0.003
text_struct/qwen3/ecfp_binary       mean      0.401     0.398  0.003
text_struct/qwen3/rdkit_descriptors mean      0.401     0.398  0.003
text_struct/f2llm/murcko_scaffold   mean      0.401     0.398  0.003
text_struct/qwen3/combined          mean      0.401     0.398  0.003
text_struct/f2llm/rdkit_descriptors mean      0.401     0.398  0.003
text_struct/bekko/rdkit_descriptors mean      0.401     0.398  0.003
text_struct/f2llm/functional_groups mean      0.401     0.398  0.003
text_struct/f2llm/iupac             mean      0.401     0.398  0.003
text_struct/f2llm/combined          mean      0.401     0.398  0.003
text_struct/qwen3/iupac             mean      0.401     0.398  0.003
text_struct/f2llm/ecfp_positional   mean      0.401     0.398  0.003
text_struct/f2llm/ecfp_binary       mean      0.401     0.398  0.003
text_struct/bekko/murcko_scaffold   mean      0.401     0.398  0.003
text_struct/bekko/functional_groups mean      0.401     0.398  0.003
text_struct/bekko/iupac             mean      0.401     0.398  0.003
string/selformer                    mean      0.401     0.398  0.003
text_semantic/f2llm/curated         mean      0.401     0.398  0.003
text_semantic/qwen3/curated         mean      0.401     0.398  0.003
text_struct/bekko/ecfp_binary       mean      0.401     0.398  0.003
text_semantic/bekko/curated         mean      0.401     0.398  0.003
text_struct/f2llm/rdkit_descriptors ridge     0.399     0.396  0.003
string/chemgpt                      gp        0.397     0.394  0.003
text_struct/qwen3/combined          knn       0.369     0.367  0.002
text_semantic/bekko/curated         knn       0.417     0.415  0.001
text_semantic/f2llm/curated         mlp       0.382     0.381  0.001
                                    knn       0.397     0.396  0.001
text_struct/f2llm/combined          mlp       0.348     0.347  0.001
text_struct/qwen3/iupac             mlp       0.342     0.341  0.000
text_semantic/bekko/curated         ridge     0.429     0.428  0.000
text_struct/qwen3/murcko_scaffold   control   0.000     0.000  0.000
text_struct/qwen3/combined          control   0.000     0.000  0.000
text_struct/f2llm/rdkit_descriptors control   0.000     0.000  0.000
text_struct/qwen3/functional_groups control   0.000     0.000  0.000
text_struct/qwen3/iupac             control   0.000     0.000  0.000
text_struct/qwen3/rdkit_descriptors control   0.000     0.000  0.000
text_struct/f2llm/ecfp_binary       control   0.000     0.000  0.000
text_struct/f2llm/murcko_scaffold   control   0.000     0.000  0.000
text_semantic/f2llm/curated         control   0.000     0.000  0.000
text_struct/bekko/combined          control   0.000     0.000  0.000
text_semantic/qwen3/curated         control   0.000     0.000  0.000
text_struct/bekko/ecfp_binary       control   0.000     0.000  0.000
text_semantic/bekko/curated         control   0.000     0.000  0.000
string/chemberta                    control   0.000     0.000  0.000
string/chemgpt                      control   0.000     0.000  0.000
string/selformer                    control   0.000     0.000  0.000
text_struct/bekko/ecfp_positional   control   0.000     0.000  0.000
sand/z                              control   0.000     0.000  0.000
ecfp_raw/counts                     control   0.000     0.000  0.000
text_struct/qwen3/ecfp_binary       control   0.000     0.000  0.000
text_struct/qwen3/ecfp_positional   control   0.000     0.000  0.000
text_struct/f2llm/ecfp_positional   control   0.000     0.000  0.000
text_struct/f2llm/combined          control   0.000     0.000  0.000
text_struct/bekko/murcko_scaffold   control   0.000     0.000  0.000
text_struct/bekko/iupac             control   0.000     0.000  0.000
text_struct/bekko/functional_groups control   0.000     0.000  0.000
text_struct/bekko/rdkit_descriptors control   0.000     0.000  0.000
text_struct/f2llm/iupac             control   0.000     0.000  0.000
text_struct/f2llm/functional_groups control   0.000     0.000  0.000
string/chemgpt                      knn       0.313     0.313 -0.000
text_struct/f2llm/iupac             mlp       0.344     0.344 -0.000
text_struct/f2llm/functional_groups mlp       0.315     0.329 -0.014
text_struct/bekko/combined          random    0.137     0.164 -0.026
text_struct/bekko/ecfp_binary       random    0.137     0.164 -0.026
string/chemberta                    random    0.137     0.164 -0.026
string/selformer                    random    0.137     0.164 -0.026
text_semantic/bekko/curated         random    0.137     0.164 -0.026
sand/z                              random    0.137     0.164 -0.026
ecfp_raw/counts                     random    0.137     0.164 -0.026
text_semantic/f2llm/curated         random    0.137     0.164 -0.026
text_semantic/qwen3/curated         random    0.137     0.164 -0.026
string/chemgpt                      random    0.137     0.164 -0.026
text_struct/bekko/ecfp_positional   random    0.137     0.164 -0.026
text_struct/bekko/iupac             random    0.137     0.164 -0.026
text_struct/bekko/rdkit_descriptors random    0.137     0.164 -0.026
text_struct/f2llm/iupac             random    0.137     0.164 -0.026
text_struct/f2llm/functional_groups random    0.137     0.164 -0.026
text_struct/f2llm/ecfp_positional   random    0.137     0.164 -0.026
text_struct/f2llm/ecfp_binary       random    0.137     0.164 -0.026
text_struct/f2llm/combined          random    0.137     0.164 -0.026
text_struct/bekko/functional_groups random    0.137     0.164 -0.026
text_struct/bekko/murcko_scaffold   random    0.137     0.164 -0.026
text_struct/f2llm/rdkit_descriptors random    0.137     0.164 -0.026
text_struct/qwen3/ecfp_binary       random    0.137     0.164 -0.026
text_struct/qwen3/combined          random    0.137     0.164 -0.026
text_struct/f2llm/murcko_scaffold   random    0.137     0.164 -0.026
text_struct/qwen3/ecfp_positional   random    0.137     0.164 -0.026
text_struct/qwen3/functional_groups random    0.137     0.164 -0.026
text_struct/qwen3/iupac             random    0.137     0.164 -0.026
text_struct/qwen3/murcko_scaffold   random    0.137     0.164 -0.026
text_struct/qwen3/rdkit_descriptors random    0.137     0.164 -0.026