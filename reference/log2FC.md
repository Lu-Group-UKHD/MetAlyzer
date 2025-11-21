# Get log2FC Data

This function returns the tibble "log2FC".

## Usage

``` r
log2FC(metalyzer_se)
```

## Arguments

- metalyzer_se:

  SummarizedExperiment

## Examples

``` r
metalyzer_se <- MetAlyzer::read_webidq(file_path = MetAlyzer::load_demodata_biocrates())
#> Input checks passed. Proceeding with reading webidq file.
#> 
#> 
#>  _____ ______   _______  _________  ________  ___           ___    ___ ________  _______   ________
#> |\   _ \  _   \|\  ___ \|\___   ___\\   __  \|\  \         |\  \  /  /|\_____  \|\  ___ \ |\   __  \
#> \ \  \\\__\ \  \ \   __/\|___ \  \_\ \  \|\  \ \  \        \ \  \/  / /\|___/  /\ \   __/|\ \  \|\  \
#>  \ \  \\|__| \  \ \  \_|/__  \ \  \ \ \   __  \ \  \        \ \    / /     /  / /\ \  \_|/_\ \   _  _\
#>   \ \  \    \ \  \ \  \_|\ \  \ \  \ \ \  \ \  \ \  \____    \/   / /     /  /_/__\ \  \_|\ \ \  \\  \| 
#>    \ \__\    \ \__\ \_______\  \ \__\ \ \__\ \__\ \_______\__/   / /     |\________\ \_______\ \__\\ _\ 
#>     \|__|     \|__|\|_______|   \|__|  \|__|\|__|\|_______|\____/ /       \|_______|\|_______|\|__|\|__|
#>                                                           \|____|/
#> 
#> 
#> Info: Reading color code "FFCBD2D7" as "#CBD2D7"
#> Info: Reading color code "FFB9DE83" as "#B9DE83"
#> Info: Reading color code "FFB9DE83" as "#B9DE83"
#> Info: Reading color code "FFA28BA3" as "#A28BA3"
#> Info: Reading color code "FFA28BA3" as "#A28BA3"
#> Info: Reading color code "FFB2D1DC" as "#B2D1DC"
#> Info: Reading color code "FF7FB2C5" as "#7FB2C5"
#> Info: Reading color code "FFB2D1DC" as "#B2D1DC"
#> Info: Reading color code "FF7FB2C5" as "#7FB2C5"
#> 
#> Measured concentration values:
#> ------------------------------
#>           0%          25%          50%          75%         100% 
#>     0.000000     0.286299     1.289381     6.308854 12522.000000 
#> 
#> NAs: 762 (3.74%)
#> 
#> 
#> Measured quantification status:
#> -------------------------------
#> Valid: 15419 (75.66%)
#> LOQ: 983 (4.82%)
#> LOD: 3978 (19.52%)
#> NAs: 0 (0%)
#> 
metalyzer_se@metadata$log2FC <- readRDS(MetAlyzer::toy_diffres())
MetAlyzer::log2FC(metalyzer_se)
#>             Metabolite                     Class        log2FC         pval
#> 1            1-Met-His        Aminoacids Related -0.4109398968 1.006211e-04
#> 2                3-IAA       Indoles Derivatives -1.5992416096 3.819741e-10
#> 3                3-IPA       Indoles Derivatives -1.2735178163 7.098227e-06
#> 4            3-Met-His        Aminoacids Related -1.1686825714 8.865990e-06
#> 5                5-AVA        Aminoacids Related  1.5654390899 8.827277e-17
#> 6                 AABA        Aminoacids Related -1.3264373525 6.082784e-10
#> 7              AbsAcid                  Hormones  0.2582790878 2.770070e-01
#> 8             AconAcid          Carboxylic Acids -1.1068126141 4.135419e-09
#> 9                 ADMA        Aminoacids Related  0.0058236329 9.361087e-01
#> 10                 Ala                Aminoacids  0.7898639261 5.391713e-10
#> 11           alpha-AAA        Aminoacids Related  1.0659645747 6.717468e-12
#> 12            Anserine        Aminoacids Related  0.5532240179 1.388262e-03
#> 13                 Arg                Aminoacids  0.1962934487 6.173235e-02
#> 14                 Asn                Aminoacids  0.3654086067 8.005909e-04
#> 15                 Asp                Aminoacids  1.5863624001 5.374475e-17
#> 16            beta-Ala           Biogenic Amines  1.3133433920 3.763443e-16
#> 17             Betaine        Aminoacids Related -0.2538509259 2.652579e-02
#> 18                  C0            Acylcarnitines -0.8292874885 6.171248e-09
#> 19                 C10            Acylcarnitines -0.5911213036 2.107198e-04
#> 20               C10:1            Acylcarnitines -0.3045171892 1.937735e-03
#> 21               C10:2            Acylcarnitines -0.3569871811 6.641629e-04
#> 22                 C12            Acylcarnitines -0.5877524133 5.072045e-05
#> 23              C12-DC            Acylcarnitines -0.0364901373 6.142900e-01
#> 24               C12:1            Acylcarnitines -0.4697382528 1.635938e-04
#> 25                 C14            Acylcarnitines  0.2202547617 1.011876e-02
#> 26               C14:1            Acylcarnitines -0.7589948666 4.764446e-05
#> 27            C14:1-OH            Acylcarnitines -0.1057445079 2.846412e-01
#> 28               C14:2            Acylcarnitines -0.2324526696 3.990545e-02
#> 29            C14:2-OH            Acylcarnitines -0.1629113706 4.108548e-02
#> 30                 C16            Acylcarnitines -0.1320052341 6.418966e-02
#> 31              C16-OH            Acylcarnitines  0.1528336825 1.890075e-01
#> 32               C16:1            Acylcarnitines -0.2099064680 9.164478e-03
#> 33            C16:1-OH            Acylcarnitines  0.3491964724 7.378438e-04
#> 34               C16:2            Acylcarnitines  0.0405936645 6.389581e-01
#> 35            C16:2-OH            Acylcarnitines  0.0529543676 4.734059e-01
#> 36                 C18            Acylcarnitines -0.2427849916 1.212235e-02
#> 37               C18:1            Acylcarnitines -1.0733079034 1.306745e-12
#> 38            C18:1-OH            Acylcarnitines  0.2794322532 2.241222e-03
#> 39               C18:2            Acylcarnitines -0.4146458890 2.303571e-04
#> 40                  C2            Acylcarnitines -0.6893485823 2.185238e-07
#> 41                  C3            Acylcarnitines  1.2262672594 1.408968e-12
#> 42       C3-DC (C4-OH)            Acylcarnitines -0.0413874445 7.258752e-01
#> 43                  C4            Acylcarnitines  1.2527345576 1.133392e-12
#> 44                C4:1            Acylcarnitines  0.0529902762 7.347141e-01
#> 45                  C5            Acylcarnitines  0.2321259453 4.776433e-02
#> 46       C5-DC (C6-OH)            Acylcarnitines -0.2875683465 2.316380e-03
#> 47             C5-M-DC            Acylcarnitines -0.0080026796 9.463176e-01
#> 48     C5-OH (C3-DC-M)            Acylcarnitines  0.2806392655 5.085923e-04
#> 49             C5:1-DC            Acylcarnitines  1.6052720804 1.655022e-07
#> 50        C6 (C4:1-DC)            Acylcarnitines -0.1949705315 2.600833e-02
#> 51                C6:1            Acylcarnitines -0.1465454949 4.118507e-02
#> 52               C7-DC            Acylcarnitines -0.4801390597 5.909189e-02
#> 53                  C8            Acylcarnitines  0.4378148423 5.927931e-02
#> 54                  C9            Acylcarnitines -1.2110318401 4.902038e-10
#> 55                  CA                Bile Acids -1.1842399997 9.256729e-04
#> 56             CE 14:0        Cholesterol Esters  0.6341370738 6.130596e-06
#> 57             CE 14:1        Cholesterol Esters  0.1642641994 1.120410e-01
#> 58             CE 15:0        Cholesterol Esters  0.4864723322 2.531029e-03
#> 59             CE 15:1        Cholesterol Esters  0.5247582290 2.165775e-04
#> 60             CE 16:0        Cholesterol Esters -0.7006321840 1.430271e-06
#> 61             CE 16:1        Cholesterol Esters -0.4213567223 6.886006e-04
#> 62             CE 17:0        Cholesterol Esters  0.5174306210 7.121490e-04
#> 63             CE 17:1        Cholesterol Esters  0.4792673992 2.781628e-04
#> 64             CE 18:0        Cholesterol Esters -0.2531478923 3.344960e-02
#> 65             CE 18:1        Cholesterol Esters -0.6787357785 4.006863e-07
#> 66             CE 18:2        Cholesterol Esters -1.5959133590 3.511108e-13
#> 67             CE 18:3        Cholesterol Esters -1.4624717130 5.002775e-15
#> 68             CE 20:1        Cholesterol Esters  0.7887607298 1.984348e-03
#> 69             CE 20:3        Cholesterol Esters -1.0220958769 1.257258e-09
#> 70             CE 20:4        Cholesterol Esters -1.3410786139 2.218286e-12
#> 71             CE 20:5        Cholesterol Esters -1.2214971157 2.444743e-08
#> 72             CE 22:1        Cholesterol Esters  0.6303593622 4.548021e-04
#> 73             CE 22:2        Cholesterol Esters  0.7267516168 3.627261e-02
#> 74             CE 22:5        Cholesterol Esters -0.0415383947 7.003344e-01
#> 75             CE 22:6        Cholesterol Esters -0.2453932291 7.813752e-02
#> 76      Cer d16:1/18:0                 Ceramides  0.4803835710 8.391435e-05
#> 77      Cer d16:1/20:0                 Ceramides  0.3155472288 7.603275e-04
#> 78      Cer d16:1/22:0                 Ceramides  0.4932967831 7.186754e-07
#> 79      Cer d16:1/23:0                 Ceramides -0.1783320070 1.102171e-01
#> 80      Cer d16:1/24:0                 Ceramides  0.3068397573 5.033963e-04
#> 81      Cer d18:0/18:0          Dihydroceramides  0.8176450649 1.563674e-08
#> 82   Cer d18:0/18:0-OH          Dihydroceramides  0.0606805163 6.190837e-01
#> 83      Cer d18:0/20:0          Dihydroceramides  0.6521621848 5.691176e-07
#> 84      Cer d18:0/22:0          Dihydroceramides  0.6731507810 1.573141e-07
#> 85      Cer d18:0/24:0          Dihydroceramides  0.3591702187 3.048785e-04
#> 86      Cer d18:0/24:1          Dihydroceramides  0.6453045587 2.489556e-08
#> 87   Cer d18:0/26:1-OH          Dihydroceramides -0.1770599150 5.268443e-01
#> 88      Cer d18:1/14:0                 Ceramides  1.1213931032 4.580785e-13
#> 89      Cer d18:1/16:0                 Ceramides  1.4505538832 1.153023e-15
#> 90      Cer d18:1/18:0                 Ceramides  1.3106367388 2.676203e-14
#> 91   Cer d18:1/18:0-OH                 Ceramides  0.3729602275 1.562144e-02
#> 92      Cer d18:1/18:1                 Ceramides  0.6204787754 4.178191e-09
#> 93      Cer d18:1/20:0                 Ceramides  1.1063495617 4.505579e-15
#> 94   Cer d18:1/20:0-OH                 Ceramides  1.1610962647 1.280097e-12
#> 95      Cer d18:1/22:0                 Ceramides  1.2934543264 7.997793e-15
#> 96      Cer d18:1/23:0                 Ceramides  0.8463941367 1.387152e-11
#> 97      Cer d18:1/24:0                 Ceramides  1.3057058797 2.837475e-15
#> 98      Cer d18:1/24:1                 Ceramides  1.2755431589 1.054158e-15
#> 99      Cer d18:1/25:0                 Ceramides  1.0092906106 1.433507e-12
#> 100     Cer d18:1/26:0                 Ceramides  1.4565736371 4.503321e-16
#> 101     Cer d18:1/26:1                 Ceramides  1.4950387636 1.073416e-15
#> 102     Cer d18:2/14:0                 Ceramides  0.7947060630 4.569129e-04
#> 103     Cer d18:2/16:0                 Ceramides  1.2958986922 7.953413e-15
#> 104     Cer d18:2/18:0                 Ceramides  1.1510563291 6.911705e-14
#> 105     Cer d18:2/18:1                 Ceramides  0.0597217661 7.913031e-01
#> 106     Cer d18:2/20:0                 Ceramides  0.8041236153 1.661160e-10
#> 107     Cer d18:2/22:0                 Ceramides  1.0543141164 1.688489e-13
#> 108     Cer d18:2/23:0                 Ceramides  0.5125382562 5.319219e-07
#> 109     Cer d18:2/24:0                 Ceramides  1.1069939434 7.382837e-15
#> 110     Cer d18:2/24:1                 Ceramides  1.0665753222 6.987620e-14
#> 111    CerP d18:1/16:0                 Ceramides  0.5144843005 3.540696e-08
#> 112            Choline      Vitamins & Cofactors  0.7084641577 2.041566e-10
#> 113                Cit        Aminoacids Related -0.1170251492 2.272597e-01
#> 114           Cortisol                  Hormones  0.0259588781 6.790344e-01
#> 115          Cortisone                  Hormones -0.7129579122 4.570531e-09
#> 116         Creatinine        Aminoacids Related -0.1926431952 1.121269e-01
#> 117                Cys                Aminoacids -0.2800566468 7.721628e-03
#> 118            Cystine        Aminoacids Related -1.5514780180 1.236493e-05
#> 119                DCA                Bile Acids -1.3641826351 4.822900e-04
#> 120       DG 14:0_14:0           Diacylglycerols -0.1649994086 8.077867e-02
#> 121       DG 14:0_18:1           Diacylglycerols  0.6047673938 6.927532e-09
#> 122       DG 14:0_18:2           Diacylglycerols -0.5310287363 2.479527e-06
#> 123       DG 14:0_20:0           Diacylglycerols -0.1867339740 1.160268e-01
#> 124       DG 14:1_18:1           Diacylglycerols  0.8039077304 2.433199e-09
#> 125       DG 14:1_20:2           Diacylglycerols  0.0265393757 8.375089e-01
#> 126       DG 16:0_16:0           Diacylglycerols  0.9490857753 1.692458e-11
#> 127       DG 16:0_16:1           Diacylglycerols  1.5711091134 5.088435e-12
#> 128       DG 16:0_18:1           Diacylglycerols  0.9002640581 9.418176e-12
#> 129       DG 16:0_18:2           Diacylglycerols -0.1292715429 2.076436e-01
#> 130       DG 16:0_20:0           Diacylglycerols -0.0167789678 8.647615e-01
#> 131       DG 16:0_20:3           Diacylglycerols  1.0391301883 2.034291e-08
#> 132       DG 16:1_18:1           Diacylglycerols  0.7335843807 9.253701e-09
#> 133       DG 16:1_18:2           Diacylglycerols -0.3316113048 5.998420e-03
#> 134       DG 16:1_20:0           Diacylglycerols  0.6453506004 1.872817e-02
#> 135       DG 17:0_18:1           Diacylglycerols  1.2288320462 1.948174e-16
#> 136       DG 18:0_20:0           Diacylglycerols  0.0469660924 7.642868e-01
#> 137       DG 18:0_20:4           Diacylglycerols  1.2339058975 2.038263e-14
#> 138       DG 18:1_18:1           Diacylglycerols  0.0376145631 5.312750e-01
#> 139       DG 18:1_18:2           Diacylglycerols -1.0591675084 3.480899e-10
#> 140       DG 18:1_18:3           Diacylglycerols -0.8326876534 3.194367e-06
#> 141       DG 18:1_18:4           Diacylglycerols  0.8699924049 1.269413e-08
#> 142       DG 18:1_20:0           Diacylglycerols -0.1354253240 4.195279e-01
#> 143       DG 18:1_20:1           Diacylglycerols  0.2152535999 9.028951e-03
#> 144       DG 18:1_20:2           Diacylglycerols  0.7668664991 4.697405e-09
#> 145       DG 18:1_20:3           Diacylglycerols  0.6540864336 7.327656e-07
#> 146       DG 18:1_20:4           Diacylglycerols  0.1855733889 9.697166e-02
#> 147       DG 18:1_22:5           Diacylglycerols -0.0127649086 9.039703e-01
#> 148       DG 18:1_22:6           Diacylglycerols  0.3536969287 3.212739e-03
#> 149       DG 18:2_18:4           Diacylglycerols -0.3964186735 1.113123e-01
#> 150       DG 18:2_20:0           Diacylglycerols -0.2505366608 1.857317e-01
#> 151       DG 18:2_20:4           Diacylglycerols -0.3894759907 5.696241e-02
#> 152       DG 18:3_18:3           Diacylglycerols  0.0814550643 7.828627e-01
#> 153       DG 18:3_20:2           Diacylglycerols -0.0088755523 9.493678e-01
#> 154       DG 22:1_22:2           Diacylglycerols  0.3313366368 1.162673e-01
#> 155     DG O-14:0_18:2           Diacylglycerols -0.4126756373 3.598647e-03
#> 156              DHEAS                  Hormones -1.5934316701 6.416126e-08
#> 157         DiCA(12:0)          Carboxylic Acids  0.0991497909 1.958267e-01
#> 158         DiCA(14:0)          Carboxylic Acids  0.0555341716 5.587242e-01
#> 159           Dopamine           Biogenic Amines  0.1302347100 1.034283e-01
#> 160            FA 14:0               Fatty Acids -0.0803718004 2.605434e-01
#> 161            FA 16:0               Fatty Acids -0.0789119227 3.405115e-01
#> 162            FA 18:0               Fatty Acids  0.0103587054 9.058831e-01
#> 163            FA 18:1               Fatty Acids -1.0235093863 1.144719e-07
#> 164            FA 18:2               Fatty Acids -1.4358021306 9.554549e-08
#> 165            FA 20:1               Fatty Acids -0.3691618303 1.277791e-02
#> 166            FA 20:2               Fatty Acids -0.3310850331 1.756243e-02
#> 167            FA 20:3               Fatty Acids -0.0747110265 7.426935e-01
#> 168    FA 20:4n-6 (AA)               Fatty Acids -0.9947906308 1.612496e-06
#> 169   FA 20:5n-3 (EPA)               Fatty Acids -0.9574525086 1.188131e-06
#> 170   FA 22:6n-3 (DHA)               Fatty Acids -1.0870433323 1.195169e-06
#> 171               GABA           Biogenic Amines  1.5987210941 1.172070e-16
#> 172                GCA                Bile Acids -1.6466036859 1.973395e-07
#> 173              GCDCA                Bile Acids -1.7376557191 5.393716e-15
#> 174               GDCA                Bile Acids -1.6895886096 9.216788e-07
#> 175               GLCA                Bile Acids -1.1232918275 2.716109e-04
#> 176              GLCAS                Bile Acids -1.5577730513 4.623608e-06
#> 177                Gln                Aminoacids -1.2852058838 3.301784e-09
#> 178                Glu                Aminoacids  1.4146985191 5.528786e-19
#> 179                Gly                Aminoacids  1.0959367263 3.064551e-14
#> 180              GUDCA                Bile Acids -1.5111853139 3.076492e-05
#> 181               HArg        Aminoacids Related -1.3214974490 1.669968e-08
#> 182               HCys        Aminoacids Related  1.2813524518 5.675884e-12
#> 183 Hex-Cer d16:1/22:0         Glycosylceramides  0.2697598126 1.583569e-02
#> 184 Hex-Cer d16:1/24:0         Glycosylceramides  0.5665614114 3.787234e-06
#> 185 Hex-Cer d18:1/14:0         Glycosylceramides  1.3441439504 2.727459e-14
#> 186 Hex-Cer d18:1/16:0         Glycosylceramides  1.5163707111 5.819861e-16
#> 187 Hex-Cer d18:1/18:0         Glycosylceramides  1.3935616774 1.294737e-13
#> 188 Hex-Cer d18:1/18:1         Glycosylceramides  0.4792213302 4.286853e-05
#> 189 Hex-Cer d18:1/20:0         Glycosylceramides  1.4001849642 1.501867e-14
#> 190 Hex-Cer d18:1/22:0         Glycosylceramides  1.3605933781 8.134925e-15
#> 191 Hex-Cer d18:1/23:0         Glycosylceramides  1.0954847668 2.706679e-13
#> 192 Hex-Cer d18:1/24:0         Glycosylceramides  1.4474629312 7.491370e-16
#> 193 Hex-Cer d18:1/24:1         Glycosylceramides  1.3066243159 2.271994e-15
#> 194 Hex-Cer d18:1/26:0         Glycosylceramides  1.5017536171 3.075301e-16
#> 195 Hex-Cer d18:1/26:1         Glycosylceramides  1.5446334012 3.422970e-16
#> 196 Hex-Cer d18:2/16:0         Glycosylceramides  1.3536545382 2.695872e-15
#> 197 Hex-Cer d18:2/18:0         Glycosylceramides  0.7931876222 1.646882e-07
#> 198 Hex-Cer d18:2/20:0         Glycosylceramides  1.0544284055 6.563889e-08
#> 199 Hex-Cer d18:2/22:0         Glycosylceramides  0.9204112408 3.166372e-11
#> 200 Hex-Cer d18:2/23:0         Glycosylceramides  0.4454116920 9.907558e-05
#> 201 Hex-Cer d18:2/24:0         Glycosylceramides  1.1948513449 5.973606e-14
#> 202 Hex2Cer d18:1/14:0         Glycosylceramides  0.3782330954 4.552590e-04
#> 203 Hex2Cer d18:1/16:0         Glycosylceramides  1.1233342198 2.278100e-13
#> 204 Hex2Cer d18:1/18:0         Glycosylceramides  1.2212794598 4.608825e-13
#> 205 Hex2Cer d18:1/20:0         Glycosylceramides  1.3727221709 7.495056e-15
#> 206 Hex2Cer d18:1/22:0         Glycosylceramides  1.4294478488 1.202805e-14
#> 207 Hex2Cer d18:1/24:0         Glycosylceramides  1.5170733113 3.620048e-16
#> 208 Hex2Cer d18:1/24:1         Glycosylceramides  1.3218869681 2.369802e-15
#> 209 Hex2Cer d18:1/26:0         Glycosylceramides  1.4857451934 1.313484e-15
#> 210 Hex2Cer d18:1/26:1         Glycosylceramides  1.5006166739 3.477251e-15
#> 211 Hex3Cer d18:1/16:0         Glycosylceramides -0.6370034042 3.268764e-07
#> 212 Hex3Cer d18:1/18:0         Glycosylceramides  0.3998780634 8.465606e-05
#> 213 Hex3Cer d18:1/20:0         Glycosylceramides  0.1508824549 3.192362e-01
#> 214 Hex3Cer d18:1/22:0         Glycosylceramides -0.5979544493 2.944279e-06
#> 215 Hex3Cer d18:1/24:1         Glycosylceramides -0.5833390309 2.162463e-06
#> 216 Hex3Cer d18:1/26:1         Glycosylceramides  0.6595308041 7.764010e-06
#> 217             Hexose                    Sugars -0.7017968730 2.130297e-04
#> 218            HipAcid          Carboxylic Acids -1.1531294494 1.350494e-05
#> 219                His                Aminoacids -0.1599309177 5.938042e-03
#> 220          Histamine           Biogenic Amines -0.4129850601 5.396767e-05
#> 221       Hypoxanthine       Nucleobases Related  0.4105160801 3.659152e-02
#> 222                Ile                Aminoacids  0.5972051710 1.035140e-09
#> 223            Ind-SO4       Indoles Derivatives -1.6074369647 6.255528e-07
#> 224             Indole       Indoles Derivatives -0.0187393942 8.210127e-01
#> 225         Kynurenine        Aminoacids Related -1.5471757805 7.258109e-14
#> 226                Lac          Carboxylic Acids -0.9372915080 2.287968e-08
#> 227                Leu                Aminoacids  0.2669943026 2.196536e-04
#> 228           LPC 14:0      Phosphatidylcholines -0.4391727047 3.271373e-06
#> 229           LPC 16:0      Phosphatidylcholines -1.2664747176 3.636851e-14
#> 230           LPC 16:1      Phosphatidylcholines -0.0710140502 2.983603e-01
#> 231           LPC 17:0      Phosphatidylcholines -0.6078349216 9.174969e-07
#> 232           LPC 18:0      Phosphatidylcholines -1.1213052802 4.111568e-13
#> 233           LPC 18:1      Phosphatidylcholines -0.7507146285 5.688671e-08
#> 234           LPC 18:2      Phosphatidylcholines -1.5734289868 9.378623e-13
#> 235           LPC 20:3      Phosphatidylcholines -0.9879050986 5.287040e-12
#> 236           LPC 20:4      Phosphatidylcholines -1.4039522487 2.795088e-13
#> 237           LPC 24:0      Phosphatidylcholines  0.5568204067 2.055654e-05
#> 238           LPC 26:0      Phosphatidylcholines  0.6833090653 1.678411e-06
#> 239           LPC 26:1      Phosphatidylcholines  0.4231715821 4.282910e-03
#> 240           LPE 12:0 Phosphatidylethanolamines -0.2358468599 7.975011e-02
#> 241           LPE 14:0 Phosphatidylethanolamines  1.2366243519 2.581176e-08
#> 242           LPE 15:0 Phosphatidylethanolamines  0.8218596767 4.502759e-09
#> 243           LPE 16:0 Phosphatidylethanolamines  0.9304817244 1.803325e-12
#> 244           LPE 16:1 Phosphatidylethanolamines  1.4076642640 1.176680e-16
#> 245           LPE 17:0 Phosphatidylethanolamines  1.2602208083 9.100286e-14
#> 246           LPE 17:1 Phosphatidylethanolamines  1.3957057070 3.078967e-11
#> 247           LPE 18:0 Phosphatidylethanolamines  1.0102541695 5.646798e-14
#> 248           LPE 18:1 Phosphatidylethanolamines  0.9949908121 2.452339e-14
#> 249           LPE 18:2 Phosphatidylethanolamines -0.6749977488 1.582287e-09
#> 250           LPE 18:3 Phosphatidylethanolamines  0.1881688677 1.005762e-01
#> 251           LPE 19:0 Phosphatidylethanolamines  0.9076817984 1.784246e-08
#> 252           LPE 19:1 Phosphatidylethanolamines  0.9274214067 2.209243e-09
#> 253           LPE 20:0 Phosphatidylethanolamines  1.1340846243 1.355300e-08
#> 254           LPE 20:1 Phosphatidylethanolamines  1.4245044280 1.217115e-11
#> 255           LPE 20:2 Phosphatidylethanolamines  1.2352551276 3.013366e-10
#> 256           LPE 20:3 Phosphatidylethanolamines  0.9630201333 2.907491e-15
#> 257           LPE 20:4 Phosphatidylethanolamines  0.4656292265 8.661607e-07
#> 258           LPE 20:5 Phosphatidylethanolamines  0.7532699764 1.234543e-05
#> 259           LPE 22:0 Phosphatidylethanolamines  0.9393862249 1.732231e-07
#> 260           LPE 22:1 Phosphatidylethanolamines  1.3528645845 2.265451e-10
#> 261           LPE 22:4 Phosphatidylethanolamines  1.1092016116 2.007831e-12
#> 262           LPE 22:5 Phosphatidylethanolamines  1.1084050437 8.041707e-14
#> 263           LPE 22:6 Phosphatidylethanolamines  0.5183809631 1.117436e-06
#> 264           LPE 24:0 Phosphatidylethanolamines  1.4128286600 3.049622e-11
#> 265         LPE P-16:0 Phosphatidylethanolamines  1.4984217356 9.075995e-17
#> 266         LPE P-18:0 Phosphatidylethanolamines  1.1684103205 3.475551e-16
#> 267         LPE P-18:1 Phosphatidylethanolamines  1.4373261649 9.731215e-17
#> 268         LPE P-20:0 Phosphatidylethanolamines  1.4120284812 3.024850e-15
#> 269         LPE P-20:4 Phosphatidylethanolamines  0.6911282606 3.492951e-07
#> 270         LPE P-22:0 Phosphatidylethanolamines  0.8604538719 1.856238e-09
#> 271           LPG 16:0     Phosphatidylglycerols  1.4267757857 6.525148e-17
#> 272           LPG 18:1     Phosphatidylglycerols  1.3516448230 4.179056e-15
#> 273           LPG 18:2     Phosphatidylglycerols  0.5855394490 1.485200e-07
#> 274           LPI 16:0     Phosphatidylinositols  1.4832036621 6.752119e-16
#> 275           LPI 18:0     Phosphatidylinositols  1.5183556054 4.505521e-18
#> 276           LPI 18:1     Phosphatidylinositols  1.4925609924 1.871820e-16
#> 277           LPI 18:2     Phosphatidylinositols -0.4267193010 4.969041e-02
#> 278           LPI 20:4     Phosphatidylinositols  0.3553054807 3.693715e-02
#> 279           LPI 22:1     Phosphatidylinositols  0.8826226344 2.713012e-12
#> 280           LPS 16:0       Phosphatidylserines  0.8237306186 5.318898e-11
#> 281           LPS 18:0       Phosphatidylserines  1.5122391619 8.488567e-17
#> 282           LPS 18:2       Phosphatidylserines  1.0932600367 1.568092e-07
#> 283           LPS 22:0       Phosphatidylserines  0.5280518599 3.312907e-05
#> 284                Lys                Aminoacids -0.2145906852 4.715582e-04
#> 285                Met                Aminoacids  0.8458655170 4.304897e-12
#> 286             Met-SO        Aminoacids Related  1.1738709521 2.776803e-07
#> 287            MG 16:1         Monoacylglycerols  0.9292850110 4.594668e-08
#> 288            MG 18:1         Monoacylglycerols  0.3371003972 4.411952e-04
#> 289            MG 18:2         Monoacylglycerols -0.2268070558 2.055295e-01
#> 290            MG 18:3         Monoacylglycerols -0.2384162820 2.392308e-01
#> 291            MG 20:1         Monoacylglycerols  0.8910376966 7.176940e-09
#> 292            MG 20:4         Monoacylglycerols  1.6124278137 5.710857e-10
#> 293            MG 20:5         Monoacylglycerols  0.0219469861 8.451781e-01
#> 294            MG 22:1         Monoacylglycerols -0.2140251085 1.403207e-01
#> 295        OH-GlutAcid          Carboxylic Acids  1.2484036137 3.386145e-16
#> 296                Orn        Aminoacids Related -1.4666913700 1.790317e-12
#> 297       p-Cresol-SO4                   Cresols -1.5836597673 1.224975e-06
#> 298       PA 17:0_18:1        Phosphatidic Acids  1.1474559344 5.283838e-13
#> 299       PA 17:1_18:1        Phosphatidic Acids  0.4664061374 5.236931e-03
#> 300       PA 18:0_18:1        Phosphatidic Acids  1.5727294768 1.805857e-16
#> 301       PA 18:0_18:2        Phosphatidic Acids  1.2629002428 4.470774e-14
#> 302       PA 18:1_18:4        Phosphatidic Acids  0.3459858281 4.476267e-02
#> 303       PA 18:1_20:0        Phosphatidic Acids  1.3305634000 1.402233e-15
#> 304       PA 18:1_20:1        Phosphatidic Acids  1.4176271697 4.426401e-16
#> 305       PA 18:1_20:2        Phosphatidic Acids  1.2523202922 8.207928e-10
#> 306       PA 18:1_20:3        Phosphatidic Acids -0.3982888460 3.072564e-02
#> 307       PA 18:1_22:0        Phosphatidic Acids  1.1285551804 4.964864e-15
#> 308       PA 18:1_22:1        Phosphatidic Acids  1.0326160550 2.035495e-06
#> 309       PA 18:1_22:3        Phosphatidic Acids -0.9617942818 3.324159e-04
#> 310       PA 18:2_20:0        Phosphatidic Acids  0.1046698588 2.054323e-01
#> 311       PA 18:2_20:1        Phosphatidic Acids  0.8161667864 1.886578e-11
#> 312       PA 18:2_22:0        Phosphatidic Acids -0.1619673572 8.397490e-02
#> 313       PA 18:2_22:1        Phosphatidic Acids -0.4889872545 3.311155e-03
#> 314       PA 20:0_20:4        Phosphatidic Acids  0.4389021743 1.115939e-07
#> 315                PAG        Aminoacids Related  1.3134674353 3.812440e-12
#> 316            PC 24:0      Phosphatidylcholines -0.1886762373 2.994892e-01
#> 317            PC 26:0      Phosphatidylcholines  0.2924960430 9.634925e-03
#> 318            PC 28:1      Phosphatidylcholines  0.1820676520 7.033153e-02
#> 319            PC 30:0      Phosphatidylcholines  1.4290799850 2.068257e-16
#> 320            PC 30:2      Phosphatidylcholines  1.4371593328 4.108978e-10
#> 321            PC 32:0      Phosphatidylcholines  1.2980967997 3.020083e-15
#> 322            PC 32:1      Phosphatidylcholines  1.4747845689 1.869057e-16
#> 323            PC 32:2      Phosphatidylcholines  1.4329934289 6.675424e-16
#> 324            PC 32:3      Phosphatidylcholines  1.1808536885 5.004436e-13
#> 325            PC 34:1      Phosphatidylcholines  0.6421069406 2.534978e-09
#> 326            PC 34:2      Phosphatidylcholines -0.2764796934 1.117470e-03
#> 327            PC 34:3      Phosphatidylcholines  0.2887492242 9.621571e-04
#> 328            PC 34:4      Phosphatidylcholines  0.3907608995 2.750142e-06
#> 329            PC 36:0      Phosphatidylcholines  0.6902587791 3.279805e-07
#> 330            PC 36:1      Phosphatidylcholines  0.6817571901 3.715423e-10
#> 331            PC 36:2      Phosphatidylcholines -0.0257915047 7.107523e-01
#> 332            PC 36:3      Phosphatidylcholines -0.6759527296 1.055497e-09
#> 333            PC 36:4      Phosphatidylcholines -1.2979791499 1.032822e-13
#> 334            PC 36:5      Phosphatidylcholines -0.8357983268 9.211218e-08
#> 335            PC 36:6      Phosphatidylcholines  0.4618584097 2.343665e-05
#> 336            PC 38:0      Phosphatidylcholines  0.3347144625 1.512452e-03
#> 337            PC 38:1      Phosphatidylcholines  0.9508026240 7.000649e-10
#> 338            PC 38:3      Phosphatidylcholines -0.5579067527 1.540521e-07
#> 339            PC 38:4      Phosphatidylcholines -1.1134445060 3.941203e-11
#> 340            PC 38:5      Phosphatidylcholines -0.8496973451 8.432561e-11
#> 341            PC 38:6      Phosphatidylcholines -1.0831949164 5.686354e-09
#> 342            PC 40:1      Phosphatidylcholines  1.0498168331 9.755325e-14
#> 343            PC 40:2      Phosphatidylcholines  1.2090106108 2.084496e-16
#> 344            PC 40:3      Phosphatidylcholines  0.7997135276 7.237239e-12
#> 345            PC 40:4      Phosphatidylcholines -0.2020905166 3.063940e-02
#> 346            PC 40:5      Phosphatidylcholines -0.3640467938 7.012655e-04
#> 347            PC 40:6      Phosphatidylcholines -0.6722996798 1.061380e-05
#> 348            PC 42:0      Phosphatidylcholines  0.1891191128 3.389164e-02
#> 349            PC 42:1      Phosphatidylcholines  1.1788854812 2.850335e-15
#> 350            PC 42:2      Phosphatidylcholines  1.2113856642 4.033110e-15
#> 351            PC 42:4      Phosphatidylcholines  0.6717943823 2.109228e-09
#> 352            PC 42:5      Phosphatidylcholines  0.4130498960 5.639636e-05
#> 353            PC 42:6      Phosphatidylcholines  0.4527573426 1.091763e-05
#> 354          PC O-28:0      Phosphatidylcholines  1.4067781042 1.730437e-15
#> 355          PC O-28:1      Phosphatidylcholines  0.7297937488 8.091244e-09
#> 356          PC O-30:0      Phosphatidylcholines  1.5348438612 8.961655e-17
#> 357          PC O-30:1      Phosphatidylcholines  1.4780287474 5.935015e-16
#> 358          PC O-30:2      Phosphatidylcholines  1.0192772906 7.176440e-12
#> 359          PC O-32:1      Phosphatidylcholines  1.4951439502 4.939022e-16
#> 360          PC O-32:2      Phosphatidylcholines  1.4047429671 1.608428e-15
#> 361          PC O-34:0      Phosphatidylcholines  1.3711609371 3.408152e-15
#> 362          PC O-34:1      Phosphatidylcholines  1.3909783776 1.528281e-15
#> 363          PC O-34:2      Phosphatidylcholines  0.9933734532 3.133678e-12
#> 364          PC O-34:3      Phosphatidylcholines -0.1914341742 4.463762e-02
#> 365          PC O-36:0      Phosphatidylcholines  1.0025972817 1.106288e-11
#> 366          PC O-36:1      Phosphatidylcholines  1.0157696591 3.183976e-14
#> 367          PC O-36:2      Phosphatidylcholines  0.8527033185 2.954852e-10
#> 368          PC O-36:3      Phosphatidylcholines  0.3870964384 1.462161e-05
#> 369          PC O-36:4      Phosphatidylcholines -0.5979606954 1.021664e-07
#> 370          PC O-36:5      Phosphatidylcholines -0.4961701808 3.074421e-06
#> 371          PC O-38:0      Phosphatidylcholines  0.3981471644 4.921028e-05
#> 372          PC O-38:2      Phosphatidylcholines  0.9902080260 9.043859e-16
#> 373          PC O-38:3      Phosphatidylcholines  0.2012889214 8.392763e-04
#> 374          PC O-38:4      Phosphatidylcholines -0.6880092640 2.300446e-09
#> 375          PC O-38:5      Phosphatidylcholines -0.1892828353 1.834000e-02
#> 376          PC O-38:6      Phosphatidylcholines  0.1053540092 2.359300e-01
#> 377          PC O-40:1      Phosphatidylcholines  0.2033535050 5.643579e-03
#> 378          PC O-40:2      Phosphatidylcholines  0.2911077662 4.140073e-04
#> 379          PC O-40:3      Phosphatidylcholines  0.3838360222 3.065543e-05
#> 380          PC O-40:4      Phosphatidylcholines -0.0996044550 9.841058e-02
#> 381          PC O-40:5      Phosphatidylcholines  0.0871547035 1.291182e-01
#> 382          PC O-40:6      Phosphatidylcholines  0.2344370754 1.912327e-02
#> 383          PC O-42:0      Phosphatidylcholines  0.4278106459 3.175312e-06
#> 384          PC O-42:1      Phosphatidylcholines  0.5011989860 5.996807e-09
#> 385          PC O-42:2      Phosphatidylcholines  0.2448970184 3.713977e-03
#> 386          PC O-42:3      Phosphatidylcholines -0.1332538782 8.024577e-02
#> 387          PC O-42:4      Phosphatidylcholines -0.0933882523 1.608324e-01
#> 388          PC O-42:5      Phosphatidylcholines  0.0109624222 8.768666e-01
#> 389          PC O-44:3      Phosphatidylcholines  0.4575817777 1.073982e-06
#> 390          PC O-44:4      Phosphatidylcholines -0.2340173705 2.039843e-03
#> 391          PC O-44:5      Phosphatidylcholines -1.0477344773 7.343384e-12
#> 392          PC O-44:6      Phosphatidylcholines -0.5297251478 7.876428e-06
#> 393            PE 20:0 Phosphatidylethanolamines  0.6927830268 6.947147e-05
#> 394            PE 28:0 Phosphatidylethanolamines  1.4873449021 7.483984e-17
#> 395            PE 28:1 Phosphatidylethanolamines  1.6151907378 8.196412e-10
#> 396            PE 30:0 Phosphatidylethanolamines  1.5475109177 9.871459e-17
#> 397            PE 30:1 Phosphatidylethanolamines  1.5500502318 1.530472e-16
#> 398            PE 31:0 Phosphatidylethanolamines  1.5706360558 1.381646e-16
#> 399            PE 32:0 Phosphatidylethanolamines  1.5372045243 1.818331e-16
#> 400            PE 32:1 Phosphatidylethanolamines  1.6028064282 1.388789e-16
#> 401            PE 32:2 Phosphatidylethanolamines  1.5959597779 2.201669e-16
#> 402            PE 33:0 Phosphatidylethanolamines  1.5541085663 1.863819e-16
#> 403            PE 33:1 Phosphatidylethanolamines  1.5929955781 4.182545e-17
#> 404            PE 33:2 Phosphatidylethanolamines  1.5701459291 9.035019e-17
#> 405            PE 34:0 Phosphatidylethanolamines  1.4439187829 9.795460e-14
#> 406            PE 34:1 Phosphatidylethanolamines  1.5853163776 2.931388e-17
#> 407            PE 34:2 Phosphatidylethanolamines  1.5121556002 2.976811e-17
#> 408            PE 34:3 Phosphatidylethanolamines  1.5563783087 1.978886e-17
#> 409            PE 34:4 Phosphatidylethanolamines  1.5800206250 5.813053e-17
#> 410            PE 35:1 Phosphatidylethanolamines  1.5906539640 2.848974e-17
#> 411            PE 35:2 Phosphatidylethanolamines  1.5511341854 2.823324e-17
#> 412            PE 35:3 Phosphatidylethanolamines  1.3857087864 3.708003e-16
#> 413            PE 36:0 Phosphatidylethanolamines  1.4455909244 4.676022e-16
#> 414            PE 36:1 Phosphatidylethanolamines  1.5337378443 2.595890e-18
#> 415            PE 36:2 Phosphatidylethanolamines  1.4224528037 1.426444e-18
#> 416            PE 36:3 Phosphatidylethanolamines  1.4344433091 6.230755e-18
#> 417            PE 36:4 Phosphatidylethanolamines  1.4025996260 6.593156e-17
#> 418            PE 36:5 Phosphatidylethanolamines  1.5059740583 2.288181e-17
#> 419            PE 36:6 Phosphatidylethanolamines  1.5492198075 1.546743e-16
#> 420            PE 38:0 Phosphatidylethanolamines  1.2450866730 1.488695e-14
#> 421            PE 38:1 Phosphatidylethanolamines  1.3284851900 3.409386e-16
#> 422            PE 38:2 Phosphatidylethanolamines  1.4560943544 4.621781e-17
#> 423            PE 38:3 Phosphatidylethanolamines  1.5170038006 1.427968e-17
#> 424            PE 38:4 Phosphatidylethanolamines  1.4013291470 2.322821e-17
#> 425            PE 38:5 Phosphatidylethanolamines  1.4627932177 2.217642e-17
#> 426            PE 38:6 Phosphatidylethanolamines  1.3128109122 2.609931e-15
#> 427            PE 38:7 Phosphatidylethanolamines  1.4776231384 3.871191e-17
#> 428            PE 40:1 Phosphatidylethanolamines  1.5284772801 6.086157e-16
#> 429            PE 40:3 Phosphatidylethanolamines  1.5128578776 1.557735e-16
#> 430            PE 40:4 Phosphatidylethanolamines  1.5468904993 2.002961e-16
#> 431            PE 40:5 Phosphatidylethanolamines  1.5384556967 6.644655e-17
#> 432            PE 40:6 Phosphatidylethanolamines  1.4729037804 1.355510e-16
#> 433            PE 40:7 Phosphatidylethanolamines  1.4899139787 9.273448e-17
#> 434            PE 40:8 Phosphatidylethanolamines  1.4131382679 3.808078e-17
#> 435            PE 42:7 Phosphatidylethanolamines  1.5614453130 2.699432e-16
#> 436            PE 42:8 Phosphatidylethanolamines  1.4850188510 3.847061e-17
#> 437           PE 44:11 Phosphatidylethanolamines  1.2887877244 5.837685e-16
#> 438           PE 44:12 Phosphatidylethanolamines  1.3461516073 4.136946e-15
#> 439            PE 44:6 Phosphatidylethanolamines  1.2716087635 1.215356e-14
#> 440            PE 44:7 Phosphatidylethanolamines  1.1540225555 7.340470e-13
#> 441     PE P-16:0/14:0 Phosphatidylethanolamines  1.4440447515 3.006695e-16
#> 442     PE P-16:0/15:0 Phosphatidylethanolamines  1.1185270156 5.812896e-16
#> 443     PE P-16:0/16:0 Phosphatidylethanolamines  1.5480368654 6.946857e-17
#> 444     PE P-16:0/16:1 Phosphatidylethanolamines  1.5912096532 2.760131e-15
#> 445     PE P-16:0/18:1 Phosphatidylethanolamines  1.5828404076 2.418718e-16
#> 446     PE P-16:0/18:2 Phosphatidylethanolamines  1.3476767790 2.638278e-14
#> 447     PE P-16:0/18:3 Phosphatidylethanolamines  1.5404121910 1.428555e-15
#> 448     PE P-16:0/20:3 Phosphatidylethanolamines  1.5529020809 1.660918e-16
#> 449     PE P-16:0/20:4 Phosphatidylethanolamines  1.4447572433 1.662350e-16
#> 450     PE P-16:0/20:5 Phosphatidylethanolamines  1.5128633742 8.821192e-17
#> 451     PE P-16:0/22:4 Phosphatidylethanolamines  1.5271048516 1.281697e-16
#> 452     PE P-16:0/22:5 Phosphatidylethanolamines  1.5191772273 3.399547e-17
#> 453     PE P-16:0/22:6 Phosphatidylethanolamines  1.4590276642 3.901327e-16
#> 454     PE P-18:0/14:0 Phosphatidylethanolamines  1.4654822679 5.997709e-17
#> 455     PE P-18:0/16:0 Phosphatidylethanolamines  1.4502241458 5.461448e-16
#> 456     PE P-18:0/16:1 Phosphatidylethanolamines  1.5551136952 2.269678e-14
#> 457     PE P-18:0/17:1 Phosphatidylethanolamines  1.5325709481 2.262828e-15
#> 458     PE P-18:0/18:0 Phosphatidylethanolamines  0.7806788940 1.748427e-04
#> 459     PE P-18:0/18:1 Phosphatidylethanolamines  1.4196867198 5.337282e-15
#> 460     PE P-18:0/18:2 Phosphatidylethanolamines  0.3838008664 7.445957e-05
#> 461     PE P-18:0/18:3 Phosphatidylethanolamines  1.2936595386 1.703452e-14
#> 462     PE P-18:0/20:1 Phosphatidylethanolamines  1.4990826518 1.843555e-15
#> 463     PE P-18:0/20:2 Phosphatidylethanolamines  1.4875220504 2.559801e-15
#> 464     PE P-18:0/20:3 Phosphatidylethanolamines  1.2775342339 2.232048e-15
#> 465     PE P-18:0/20:4 Phosphatidylethanolamines  0.7409088929 2.285894e-12
#> 466     PE P-18:0/20:5 Phosphatidylethanolamines  0.9407414684 2.131814e-10
#> 467     PE P-18:0/22:3 Phosphatidylethanolamines  1.4993914255 8.841907e-18
#> 468     PE P-18:0/22:4 Phosphatidylethanolamines  1.2706676162 1.470682e-17
#> 469     PE P-18:0/22:5 Phosphatidylethanolamines  1.3081400154 2.500452e-17
#> 470     PE P-18:0/22:6 Phosphatidylethanolamines  1.0534696078 6.469857e-13
#> 471     PE P-18:1/18:1 Phosphatidylethanolamines  1.5291813141 2.021705e-16
#> 472     PE P-18:1/18:2 Phosphatidylethanolamines  1.0399176609 6.691167e-13
#> 473     PE P-18:1/20:4 Phosphatidylethanolamines  1.1976636408 1.312011e-14
#> 474     PE P-18:1/20:5 Phosphatidylethanolamines  1.3302977097 3.948081e-15
#> 475     PE P-18:1/22:6 Phosphatidylethanolamines  1.3623139372 8.481157e-15
#> 476     PE P-20:0/14:0 Phosphatidylethanolamines  0.1394308946 1.808087e-01
#> 477     PE P-20:0/16:0 Phosphatidylethanolamines  1.4915961566 9.376290e-11
#> 478     PE P-20:0/16:1 Phosphatidylethanolamines  1.2721131567 4.891295e-13
#> 479     PE P-20:0/17:1 Phosphatidylethanolamines  0.6341808446 3.245065e-05
#> 480     PE P-20:0/18:1 Phosphatidylethanolamines  0.9585505859 5.249430e-12
#> 481     PE P-20:0/18:2 Phosphatidylethanolamines -0.3858586769 8.144321e-05
#> 482     PE P-20:0/20:0 Phosphatidylethanolamines  0.3020336617 1.413258e-01
#> 483     PE P-20:0/20:4 Phosphatidylethanolamines -0.0583913154 3.789501e-01
#> 484     PE P-20:0/20:5 Phosphatidylethanolamines  0.8859410037 4.784070e-09
#> 485       PG 14:0_16:0     Phosphatidylglycerols  1.3735087069 1.574533e-15
#> 486       PG 16:0_16:0     Phosphatidylglycerols  1.5987783668 4.590022e-12
#> 487       PG 16:0_16:1     Phosphatidylglycerols  1.5055890191 4.850348e-16
#> 488       PG 16:0_18:1     Phosphatidylglycerols  1.6011636234 2.688050e-16
#> 489       PG 16:0_18:2     Phosphatidylglycerols  1.0456365369 1.598515e-14
#> 490       PG 16:0_18:3     Phosphatidylglycerols -0.0008422691 9.902368e-01
#> 491       PG 16:0_19:1     Phosphatidylglycerols  1.5453174137 6.191444e-17
#> 492       PG 16:0_20:3     Phosphatidylglycerols  1.2557395395 1.453738e-15
#> 493       PG 16:0_20:4     Phosphatidylglycerols  1.2479594791 3.182192e-15
#> 494       PG 16:0_20:5     Phosphatidylglycerols  1.2938260355 2.423378e-15
#> 495       PG 16:0_22:1     Phosphatidylglycerols  0.2540404731 1.680289e-03
#> 496       PG 16:0_22:2     Phosphatidylglycerols  1.0928305457 2.574789e-14
#> 497       PG 16:1_16:1     Phosphatidylglycerols  1.6468217436 2.878218e-12
#> 498       PG 16:1_18:0     Phosphatidylglycerols  1.4441246309 5.773488e-16
#> 499       PG 16:1_18:1     Phosphatidylglycerols  1.3836953768 1.127332e-16
#> 500       PG 16:1_18:2     Phosphatidylglycerols  0.2296186165 1.984439e-04
#> 501       PG 16:1_20:4     Phosphatidylglycerols  1.0590632220 4.492445e-14
#> 502       PG 16:2_18:1     Phosphatidylglycerols  1.5250487107 3.079367e-18
#> 503       PG 16:2_18:2     Phosphatidylglycerols  1.1547572592 5.439261e-16
#> 504       PG 16:3_18:1     Phosphatidylglycerols  0.9644554181 3.423242e-13
#> 505       PG 17:1_18:1     Phosphatidylglycerols  1.5138402797 7.692109e-17
#> 506       PG 18:0_18:1     Phosphatidylglycerols  1.5441656038 7.904231e-15
#> 507       PG 18:0_18:2     Phosphatidylglycerols  1.3757527660 1.866776e-11
#> 508       PG 18:0_18:3     Phosphatidylglycerols -0.5801943341 4.202868e-08
#> 509       PG 18:0_22:1     Phosphatidylglycerols -0.2385709999 5.708648e-03
#> 510       PG 18:1_18:1     Phosphatidylglycerols  1.5590488861 1.499132e-16
#> 511       PG 18:1_18:2     Phosphatidylglycerols  0.8271527732 9.684075e-11
#> 512       PG 18:1_18:3     Phosphatidylglycerols  1.0323574194 1.105905e-13
#> 513       PG 18:1_20:0     Phosphatidylglycerols  1.1191249949 4.687125e-14
#> 514       PG 18:1_20:1     Phosphatidylglycerols  1.5276401863 9.866771e-12
#> 515       PG 18:1_20:3     Phosphatidylglycerols  0.6681116037 3.616727e-10
#> 516       PG 18:1_20:4     Phosphatidylglycerols  0.5328363892 3.912061e-08
#> 517       PG 18:1_20:5     Phosphatidylglycerols  1.3399734349 2.128530e-15
#> 518       PG 18:1_22:0     Phosphatidylglycerols  1.2738404624 8.758333e-15
#> 519       PG 18:1_22:1     Phosphatidylglycerols  1.3830252964 1.703958e-15
#> 520       PG 18:1_22:2     Phosphatidylglycerols  1.4073621683 5.088860e-18
#> 521       PG 18:1_22:3     Phosphatidylglycerols  1.0619837681 1.092550e-14
#> 522       PG 18:1_22:4     Phosphatidylglycerols  0.7730968946 2.569601e-11
#> 523       PG 18:1_22:5     Phosphatidylglycerols  1.0064885087 1.601101e-13
#> 524       PG 18:2_18:2     Phosphatidylglycerols  0.5061094024 9.550625e-07
#> 525       PG 18:2_18:3     Phosphatidylglycerols -0.3848894025 4.945425e-06
#> 526       PG 18:2_18:4     Phosphatidylglycerols -0.1535378342 8.518001e-02
#> 527       PG 18:2_20:0     Phosphatidylglycerols -0.5257454039 2.415181e-06
#> 528       PG 18:2_20:3     Phosphatidylglycerols -1.1526821106 1.998039e-13
#> 529       PG 18:2_20:4     Phosphatidylglycerols -1.2718338273 5.887546e-14
#> 530       PG 18:2_20:5     Phosphatidylglycerols  0.0827655241 2.946268e-01
#> 531       PG 18:2_22:0     Phosphatidylglycerols -0.3881605057 1.905355e-04
#> 532       PG 18:2_22:1     Phosphatidylglycerols -0.3045452160 1.134576e-03
#> 533       PG 18:2_22:3     Phosphatidylglycerols -0.9803338059 5.536389e-12
#> 534       PG 18:2_22:4     Phosphatidylglycerols -1.3227535506 9.130403e-14
#> 535       PG 20:3_20:4     Phosphatidylglycerols -0.9658042994 8.070395e-11
#> 536       PG 20:4_20:4     Phosphatidylglycerols -1.3935120286 1.121269e-13
#> 537       PG 20:4_22:1     Phosphatidylglycerols -0.1698772326 8.372552e-02
#> 538       PG 20:4_22:3     Phosphatidylglycerols -0.3895607068 3.130856e-04
#> 539       PG 20:4_22:4     Phosphatidylglycerols -1.3224507271 5.161162e-12
#> 540       PG 22:4_22:6     Phosphatidylglycerols -0.5798040149 1.169292e-05
#> 541       PG 22:5_22:6     Phosphatidylglycerols -0.7719277358 2.900625e-09
#> 542       PG 22:6_22:6     Phosphatidylglycerols -0.7933414119 1.651269e-06
#> 543                Phe                Aminoacids  0.3505162619 1.284937e-05
#> 544      PheAlaBetaine        Aminoacids Related -0.3218977669 2.246758e-01
#> 545       PI 14:0_18:1     Phosphatidylinositols  1.3624470387 5.954077e-16
#> 546       PI 14:0_18:2     Phosphatidylinositols -0.1697652164 2.403803e-02
#> 547       PI 15:0_16:0     Phosphatidylinositols  0.4451853750 9.389502e-07
#> 548       PI 15:1_16:0     Phosphatidylinositols  0.7965994580 2.139262e-11
#> 549       PI 16:0_16:0     Phosphatidylinositols  1.3725536157 2.058824e-14
#> 550       PI 16:0_17:0     Phosphatidylinositols  0.8769938149 1.193328e-09
#> 551       PI 16:0_17:1     Phosphatidylinositols  1.2698717307 1.273277e-16
#> 552       PI 16:0_17:2     Phosphatidylinositols  0.5255649520 1.166053e-07
#> 553       PI 16:0_18:1     Phosphatidylinositols  1.4241778151 6.982236e-17
#> 554       PI 16:0_18:2     Phosphatidylinositols  0.9177217112 6.739375e-11
#> 555       PI 16:0_18:3     Phosphatidylinositols  1.1827473070 7.409156e-16
#> 556       PI 16:0_20:0     Phosphatidylinositols -0.5791181158 3.499524e-05
#> 557       PI 16:0_20:3     Phosphatidylinositols  1.2230880954 2.340921e-15
#> 558       PI 16:0_20:4     Phosphatidylinositols  1.1245718656 3.950187e-13
#> 559       PI 16:0_22:1     Phosphatidylinositols  0.5536880108 3.440517e-08
#> 560       PI 16:1_18:0     Phosphatidylinositols  1.3826418092 3.437306e-15
#> 561       PI 16:1_18:1     Phosphatidylinositols  1.4358790547 2.716119e-16
#> 562       PI 16:1_18:2     Phosphatidylinositols -0.3379810866 1.229956e-04
#> 563       PI 17:1_18:1     Phosphatidylinositols  0.7479372629 1.053750e-10
#> 564       PI 17:1_18:2     Phosphatidylinositols -1.3599330684 4.234521e-14
#> 565       PI 18:0_18:0     Phosphatidylinositols  0.8687232781 4.645506e-11
#> 566       PI 18:0_18:1     Phosphatidylinositols  1.2873324421 2.488615e-15
#> 567       PI 18:0_18:2     Phosphatidylinositols  0.6447092736 1.586533e-07
#> 568       PI 18:0_18:3     Phosphatidylinositols  1.3251712858 1.817380e-14
#> 569       PI 18:0_20:0     Phosphatidylinositols -0.1966228724 8.642401e-02
#> 570       PI 18:0_20:3     Phosphatidylinositols  0.8943429160 2.228155e-12
#> 571       PI 18:0_20:4     Phosphatidylinositols  0.5725545459 2.374655e-08
#> 572       PI 18:0_22:0     Phosphatidylinositols  0.8061069017 2.870223e-10
#> 573       PI 18:1_18:1     Phosphatidylinositols  1.5278863756 2.298962e-16
#> 574       PI 18:1_18:2     Phosphatidylinositols  1.3720395059 8.397100e-16
#> 575       PI 18:1_18:3     Phosphatidylinositols  1.5441129748 4.361026e-12
#> 576       PI 18:1_20:0     Phosphatidylinositols  0.9790228254 5.067131e-12
#> 577       PI 18:1_20:1     Phosphatidylinositols  1.3472706679 1.879807e-15
#> 578       PI 18:1_20:2     Phosphatidylinositols  1.4275313462 6.433934e-16
#> 579       PI 18:1_20:3     Phosphatidylinositols  1.5813839351 1.164449e-12
#> 580       PI 18:1_20:4     Phosphatidylinositols  1.4724923151 4.815710e-16
#> 581       PI 18:1_20:5     Phosphatidylinositols  1.5766553094 7.795675e-12
#> 582       PI 18:1_22:0     Phosphatidylinositols  1.3999594018 4.898960e-16
#> 583       PI 18:1_22:1     Phosphatidylinositols  1.0174807187 5.299508e-13
#> 584       PI 18:1_22:2     Phosphatidylinositols  1.2463307460 1.223507e-14
#> 585       PI 18:1_22:3     Phosphatidylinositols  1.2857544590 7.942128e-10
#> 586       PI 18:1_22:5     Phosphatidylinositols  1.4820602373 1.202966e-11
#> 587       PI 18:1_22:6     Phosphatidylinositols  1.3901370135 6.609473e-17
#> 588       PI 18:2_20:0     Phosphatidylinositols  0.1298459970 1.509258e-01
#> 589       PI 18:2_20:1     Phosphatidylinositols  0.2587784962 1.030406e-03
#> 590       PI 18:2_20:5     Phosphatidylinositols  0.6404813591 1.518730e-09
#> 591       PI 18:2_22:0     Phosphatidylinositols  0.6546426833 9.591845e-08
#> 592       PI 18:2_22:1     Phosphatidylinositols -0.4240131436 1.096611e-04
#> 593       PI 18:2_22:6     Phosphatidylinositols  0.4422781270 3.068291e-07
#> 594                Pro                Aminoacids  0.2666480145 6.923769e-04
#> 595         ProBetaine        Aminoacids Related -0.9667113010 6.636859e-03
#> 596            PS 30:0       Phosphatidylserines  1.3774409725 2.324396e-14
#> 597            PS 34:1       Phosphatidylserines  1.6169030486 4.952466e-17
#> 598            PS 34:2       Phosphatidylserines  1.6091010976 8.388718e-17
#> 599            PS 36:1       Phosphatidylserines  1.6180994966 2.115794e-17
#> 600            PS 36:2       Phosphatidylserines  1.6202326610 3.283947e-17
#> 601            PS 36:3       Phosphatidylserines  1.6166476906 6.369295e-17
#> 602            PS 36:4       Phosphatidylserines  1.4674559867 8.064654e-16
#> 603            PS 36:5       Phosphatidylserines  0.3477995229 3.910673e-04
#> 604            PS 38:4       Phosphatidylserines  1.5642195185 6.778065e-18
#> 605            PS 38:5       Phosphatidylserines  1.4967468451 1.834913e-15
#> 606            PS 38:6       Phosphatidylserines  1.5403267332 1.555128e-15
#> 607            PS 38:7       Phosphatidylserines  0.6046346776 1.388525e-07
#> 608            PS 40:4       Phosphatidylserines  1.6150328372 6.396467e-15
#> 609            PS 40:5       Phosphatidylserines  1.6300937657 1.327878e-15
#> 610            PS 40:6       Phosphatidylserines  1.6304691471 1.670880e-15
#> 611            PS 40:7       Phosphatidylserines  1.3442782046 1.618373e-14
#> 612            PS 40:8       Phosphatidylserines  0.9437634728 5.819517e-12
#> 613         Putrescine           Biogenic Amines  0.3924580261 3.292833e-05
#> 614          Sarcosine        Aminoacids Related -1.0150231036 2.895060e-08
#> 615               SDMA        Aminoacids Related -1.0609708232 7.411151e-13
#> 616                Ser                Aminoacids  0.8452131458 2.256788e-11
#> 617          Serotonin           Biogenic Amines -1.0713877800 4.038304e-08
#> 618            SM 33:1            Sphingomyelins  0.4493265849 1.726025e-04
#> 619            SM 34:1            Sphingomyelins  0.7503922825 7.180858e-10
#> 620            SM 34:2            Sphingomyelins  0.2833312174 3.040951e-03
#> 621            SM 35:1            Sphingomyelins  0.5914312457 3.372046e-06
#> 622            SM 36:1            Sphingomyelins -0.0867388318 3.203262e-01
#> 623            SM 36:2            Sphingomyelins -0.6084721883 2.041220e-06
#> 624            SM 41:1            Sphingomyelins -0.1966484333 5.473107e-02
#> 625            SM 41:2            Sphingomyelins -0.0385838290 7.219385e-01
#> 626            SM 42:1            Sphingomyelins  0.6869972223 1.279242e-09
#> 627            SM 42:2            Sphingomyelins  0.3208966064 3.049597e-04
#> 628            SM 43:1            Sphingomyelins  0.2538856872 1.427325e-02
#> 629            SM 44:1            Sphingomyelins  1.1084259245 1.248819e-09
#> 630            SM 44:2            Sphingomyelins  1.0902752652 3.552798e-12
#> 631         SPBP d14:1           Sphingoid Bases -1.6126824850 2.032645e-09
#> 632         SPBP d16:1           Sphingoid Bases  0.0904167766 5.771551e-01
#> 633         SPBP d17:0           Sphingoid Bases  0.1035228091 3.154568e-01
#> 634         SPBP d17:1           Sphingoid Bases -1.3179063453 7.298336e-07
#> 635         Spermidine           Biogenic Amines  1.5394469462 1.002698e-15
#> 636                Suc          Carboxylic Acids  0.8786835089 2.638751e-14
#> 637          t4-OH-Pro        Aminoacids Related  0.5481790446 1.338927e-06
#> 638            Taurine        Aminoacids Related  1.1937914355 1.694946e-16
#> 639                TCA                Bile Acids -0.8445976429 7.826660e-04
#> 640              TCDCA                Bile Acids -0.8559414592 2.628129e-04
#> 641       TG 14:0_32:2          Triacylglycerols -1.2436774481 1.891718e-07
#> 642       TG 14:0_34:0          Triacylglycerols -0.8638230158 2.084631e-04
#> 643       TG 14:0_34:1          Triacylglycerols -1.1646091931 5.917695e-07
#> 644       TG 14:0_34:2          Triacylglycerols -1.4340181844 1.002581e-09
#> 645       TG 14:0_34:3          Triacylglycerols -1.6711127412 1.889543e-09
#> 646       TG 14:0_35:1          Triacylglycerols -0.5377592957 4.991344e-03
#> 647       TG 14:0_35:2          Triacylglycerols -0.8629863820 2.708372e-06
#> 648       TG 14:0_36:1          Triacylglycerols -0.9799342969 8.107648e-05
#> 649       TG 14:0_36:2          Triacylglycerols -1.5851491612 4.026517e-12
#> 650       TG 14:0_36:3          Triacylglycerols -1.6897734658 1.139206e-11
#> 651       TG 14:0_38:4          Triacylglycerols -1.1445110987 3.583239e-09
#> 652       TG 14:0_38:5          Triacylglycerols -1.2287158280 3.853742e-08
#> 653       TG 14:0_39:3          Triacylglycerols -0.3245087934 1.562991e-01
#> 654       TG 16:0_28:1          Triacylglycerols -1.2500579434 4.550822e-05
#> 655       TG 16:0_28:2          Triacylglycerols -1.3727451790 4.319476e-05
#> 656       TG 16:0_30:2          Triacylglycerols -1.5052546332 1.899379e-07
#> 657       TG 16:0_32:0          Triacylglycerols -0.2591642821 1.895767e-01
#> 658       TG 16:0_32:1          Triacylglycerols -0.7928966128 6.149997e-05
#> 659       TG 16:0_32:2          Triacylglycerols -1.2924388980 3.405146e-08
#> 660       TG 16:0_32:3          Triacylglycerols -1.5339813270 1.753824e-08
#> 661       TG 16:0_33:1          Triacylglycerols -0.2538995209 9.409450e-02
#> 662       TG 16:0_33:2          Triacylglycerols -0.7994561230 9.026690e-06
#> 663       TG 16:0_34:0          Triacylglycerols -0.6612907231 7.075455e-03
#> 664       TG 16:0_34:1          Triacylglycerols -1.2154823523 5.584381e-07
#> 665       TG 16:0_34:2          Triacylglycerols -1.4968803078 4.452336e-10
#> 666       TG 16:0_34:3          Triacylglycerols -1.7327937476 4.681240e-10
#> 667       TG 16:0_34:4          Triacylglycerols -1.6641701218 1.055756e-08
#> 668       TG 16:0_35:1          Triacylglycerols -0.8100536062 1.088861e-05
#> 669       TG 16:0_35:2          Triacylglycerols -1.1579290847 2.506144e-10
#> 670       TG 16:0_35:3          Triacylglycerols -1.5894974367 4.232577e-09
#> 671       TG 16:0_36:2          Triacylglycerols -1.7749395265 4.185211e-11
#> 672       TG 16:0_36:4          Triacylglycerols -1.7143206828 1.525836e-10
#> 673       TG 16:0_36:6          Triacylglycerols -1.5034598951 1.234287e-07
#> 674       TG 16:0_37:3          Triacylglycerols -1.2238830658 1.046573e-06
#> 675       TG 16:0_38:1          Triacylglycerols -0.8971714854 3.861653e-05
#> 676       TG 16:0_38:2          Triacylglycerols -1.2666124120 3.510674e-10
#> 677       TG 16:0_38:3          Triacylglycerols -1.4340729754 3.826442e-12
#> 678       TG 16:0_38:4          Triacylglycerols -1.4902493735 4.616011e-12
#> 679       TG 16:0_38:5          Triacylglycerols -1.5260917859 3.177656e-11
#> 680       TG 16:0_38:6          Triacylglycerols -1.5180604118 2.845603e-11
#> 681       TG 16:0_38:7          Triacylglycerols -1.5400631410 5.562700e-09
#> 682       TG 16:0_40:6          Triacylglycerols -1.6989217546 6.905062e-12
#> 683       TG 16:0_40:7          Triacylglycerols -1.6233350000 1.766420e-14
#> 684       TG 16:1_28:0          Triacylglycerols -0.6148956574 2.590616e-03
#> 685       TG 16:1_30:1          Triacylglycerols -0.7691267335 2.054073e-04
#> 686       TG 16:1_32:0          Triacylglycerols -0.1958406336 2.569213e-01
#> 687       TG 16:1_32:1          Triacylglycerols -0.7625582743 9.722187e-05
#> 688       TG 16:1_32:2          Triacylglycerols -1.1946819133 1.721557e-08
#> 689       TG 16:1_33:1          Triacylglycerols -0.1921904569 1.452541e-01
#> 690       TG 16:1_34:0          Triacylglycerols -0.3659163756 5.040693e-02
#> 691       TG 16:1_34:1          Triacylglycerols -1.2826364156 9.079636e-09
#> 692       TG 16:1_34:2          Triacylglycerols -1.4777894390 3.448148e-10
#> 693       TG 16:1_34:3          Triacylglycerols -1.6877300391 1.830055e-09
#> 694       TG 16:1_36:1          Triacylglycerols -0.9321303088 1.055620e-06
#> 695       TG 16:1_36:2          Triacylglycerols -1.5604716604 2.749500e-12
#> 696       TG 16:1_38:3          Triacylglycerols -1.2413002725 1.161141e-09
#> 697       TG 16:1_38:4          Triacylglycerols -1.3751369253 9.658585e-11
#> 698       TG 16:1_38:5          Triacylglycerols -1.4181968661 3.369172e-10
#> 699       TG 17:0_32:1          Triacylglycerols -0.3358234457 2.317392e-02
#> 700       TG 17:0_34:2          Triacylglycerols -1.4093254220 5.443162e-09
#> 701       TG 17:1_32:1          Triacylglycerols -0.6174139768 1.346812e-04
#> 702       TG 17:1_34:1          Triacylglycerols -1.1767097059 2.964426e-09
#> 703       TG 17:1_34:2          Triacylglycerols -1.4003335021 4.049864e-12
#> 704       TG 17:1_34:3          Triacylglycerols -1.5596667881 1.958832e-09
#> 705       TG 17:1_36:3          Triacylglycerols -1.6258935900 9.076004e-13
#> 706       TG 17:1_36:5          Triacylglycerols -1.1309348067 9.842510e-06
#> 707       TG 17:1_38:5          Triacylglycerols -0.5257339627 3.830009e-04
#> 708       TG 17:1_38:6          Triacylglycerols -0.4541915839 4.211565e-03
#> 709       TG 17:2_34:2          Triacylglycerols -1.0838317696 3.209962e-06
#> 710       TG 17:2_34:3          Triacylglycerols -1.3010067257 1.255656e-06
#> 711       TG 17:2_36:2          Triacylglycerols -0.5363360490 3.209595e-05
#> 712       TG 17:2_36:3          Triacylglycerols -1.2773333228 7.593981e-08
#> 713       TG 17:2_36:4          Triacylglycerols -1.0557120611 9.561649e-06
#> 714       TG 17:2_38:5          Triacylglycerols -0.2736696610 5.746228e-02
#> 715       TG 17:2_38:6          Triacylglycerols -0.2451227740 2.968135e-01
#> 716       TG 18:0_30:0          Triacylglycerols -0.5897724246 1.631983e-02
#> 717       TG 18:0_30:1          Triacylglycerols -0.9307052885 2.593604e-04
#> 718       TG 18:0_32:0          Triacylglycerols -0.6216336709 1.171921e-02
#> 719       TG 18:0_32:1          Triacylglycerols -0.7748527311 4.304380e-04
#> 720       TG 18:0_32:2          Triacylglycerols -1.1508886354 2.522312e-06
#> 721       TG 18:0_34:2          Triacylglycerols -1.3764620499 2.789261e-08
#> 722       TG 18:0_34:3          Triacylglycerols -1.5623507002 4.360575e-09
#> 723       TG 18:0_36:1          Triacylglycerols -0.9772784925 1.793205e-03
#> 724       TG 18:0_36:2          Triacylglycerols -1.4293031939 2.204974e-07
#> 725       TG 18:0_36:3          Triacylglycerols -1.6387328777 1.671995e-11
#> 726       TG 18:0_36:4          Triacylglycerols -1.5952482127 3.829149e-10
#> 727       TG 18:0_36:5          Triacylglycerols -1.5856553580 6.920518e-09
#> 728       TG 18:0_38:6          Triacylglycerols -1.5857172164 6.675251e-11
#> 729       TG 18:1_26:0          Triacylglycerols -1.4236185016 1.724588e-05
#> 730       TG 18:1_28:1          Triacylglycerols -1.4701144664 4.542087e-07
#> 731       TG 18:1_30:0          Triacylglycerols -1.2386765067 3.656900e-07
#> 732       TG 18:1_30:1          Triacylglycerols -1.4300554951 6.219814e-10
#> 733       TG 18:1_30:2          Triacylglycerols -1.5698686712 5.280209e-10
#> 734       TG 18:1_31:0          Triacylglycerols -0.5751229197 7.183069e-04
#> 735       TG 18:1_32:0          Triacylglycerols -1.2479296856 1.913381e-07
#> 736       TG 18:1_32:1          Triacylglycerols -1.5453987222 1.566729e-09
#> 737       TG 18:1_32:2          Triacylglycerols -1.6027198878 8.086557e-13
#> 738       TG 18:1_32:3          Triacylglycerols -1.6446519941 3.813561e-11
#> 739       TG 18:1_33:0          Triacylglycerols -0.9961155233 4.422028e-07
#> 740       TG 18:1_33:1          Triacylglycerols -1.2112202730 7.448001e-11
#> 741       TG 18:1_33:2          Triacylglycerols -1.3624272039 6.655116e-12
#> 742       TG 18:1_33:3          Triacylglycerols -1.4408187277 4.029657e-09
#> 743       TG 18:1_34:1          Triacylglycerols -1.6443053261 9.068672e-13
#> 744       TG 18:1_34:2          Triacylglycerols -1.7081227879 1.481712e-14
#> 745       TG 18:1_34:3          Triacylglycerols -1.7120970131 7.927328e-14
#> 746       TG 18:1_35:2          Triacylglycerols -1.6134400771 2.017801e-11
#> 747       TG 18:1_35:3          Triacylglycerols -1.6168173483 1.441991e-10
#> 748       TG 18:1_36:0          Triacylglycerols -1.3066240057 2.281488e-05
#> 749       TG 18:1_36:1          Triacylglycerols -1.5065461919 3.724109e-07
#> 750       TG 18:1_36:2          Triacylglycerols -1.6885425376 4.890184e-12
#> 751       TG 18:1_36:4          Triacylglycerols -1.7659470269 1.534018e-08
#> 752       TG 18:1_38:5          Triacylglycerols -1.6213382364 1.598295e-15
#> 753       TG 18:1_38:6          Triacylglycerols -1.5980310434 1.621146e-14
#> 754       TG 18:2_30:0          Triacylglycerols -1.6745762109 9.397100e-08
#> 755       TG 18:2_31:0          Triacylglycerols -0.8107448558 1.296528e-07
#> 756       TG 18:2_32:1          Triacylglycerols -1.6710441756 2.903261e-12
#> 757       TG 18:2_32:2          Triacylglycerols -1.6866010300 1.025275e-10
#> 758       TG 18:2_33:0          Triacylglycerols -1.5855313682 2.256852e-07
#> 759       TG 18:2_33:1          Triacylglycerols -1.6731752473 1.005721e-09
#> 760       TG 18:2_33:2          Triacylglycerols -1.5421742340 5.735638e-09
#> 761       TG 18:2_34:0          Triacylglycerols -1.5826946086 1.487128e-07
#> 762       TG 18:2_34:1          Triacylglycerols -1.7161859666 6.166288e-13
#> 763       TG 18:2_34:2          Triacylglycerols -1.8090089054 6.600468e-10
#> 764       TG 18:2_34:4          Triacylglycerols -1.7399764012 9.626123e-10
#> 765       TG 18:2_35:2          Triacylglycerols -1.7301270528 1.287798e-09
#> 766       TG 18:2_35:3          Triacylglycerols -1.6271829587 1.733238e-08
#> 767       TG 18:2_36:0          Triacylglycerols -1.3787476882 4.163865e-07
#> 768       TG 18:2_36:1          Triacylglycerols -1.7583994951 3.765261e-10
#> 769       TG 18:2_38:5          Triacylglycerols -1.7793651905 1.299790e-10
#> 770       TG 18:2_38:6          Triacylglycerols -1.5920208216 4.972298e-11
#> 771       TG 18:3_30:0          Triacylglycerols -1.5051437698 1.607045e-07
#> 772       TG 18:3_32:0          Triacylglycerols -1.5341911936 5.509980e-08
#> 773       TG 18:3_32:1          Triacylglycerols -1.6137981450 4.446379e-10
#> 774       TG 18:3_34:0          Triacylglycerols -1.6232765923 3.030904e-08
#> 775       TG 18:3_38:6          Triacylglycerols -1.5539968678 4.886746e-09
#> 776       TG 20:1_24:3          Triacylglycerols -0.8031331449 3.659395e-06
#> 777       TG 20:1_26:1          Triacylglycerols -0.7678761428 1.711742e-04
#> 778       TG 20:1_30:1          Triacylglycerols -0.9330745743 1.635537e-05
#> 779       TG 20:1_31:0          Triacylglycerols -0.3854058505 7.127523e-04
#> 780       TG 20:1_32:1          Triacylglycerols -0.8324914281 3.274667e-06
#> 781       TG 20:1_32:2          Triacylglycerols -1.2955163850 3.158425e-07
#> 782       TG 20:1_32:3          Triacylglycerols -1.2959185077 7.706453e-10
#> 783       TG 20:1_34:1          Triacylglycerols -1.3572327708 2.678629e-10
#> 784       TG 20:1_34:2          Triacylglycerols -1.5100101133 3.349416e-11
#> 785       TG 20:2_32:0          Triacylglycerols -0.8679152819 2.068809e-05
#> 786       TG 20:2_32:1          Triacylglycerols -1.3212988959 2.031187e-10
#> 787       TG 20:2_34:1          Triacylglycerols -1.5258838350 1.030637e-13
#> 788       TG 20:2_34:2          Triacylglycerols -1.7326126986 7.486547e-11
#> 789       TG 20:2_34:4          Triacylglycerols -0.9621070426 7.961505e-06
#> 790       TG 20:3_32:0          Triacylglycerols -1.0793919811 4.495384e-08
#> 791       TG 20:3_32:1          Triacylglycerols -1.3963446368 1.480708e-10
#> 792       TG 20:3_32:2          Triacylglycerols -1.5876384177 1.839921e-09
#> 793       TG 20:3_34:0          Triacylglycerols -0.1891431964 3.441201e-01
#> 794       TG 20:3_34:1          Triacylglycerols -1.6889881787 8.833757e-11
#> 795       TG 20:3_34:2          Triacylglycerols -1.7648460795 1.051415e-11
#> 796       TG 20:3_34:3          Triacylglycerols -1.5888253726 3.008789e-12
#> 797       TG 20:3_36:4          Triacylglycerols -1.4032398822 3.000732e-10
#> 798       TG 20:3_36:5          Triacylglycerols -1.3684532638 6.112143e-09
#> 799       TG 20:4_30:0          Triacylglycerols -1.3698165387 1.655740e-06
#> 800       TG 20:4_32:0          Triacylglycerols -1.2892251153 1.634385e-06
#> 801       TG 20:4_32:1          Triacylglycerols -1.4370211951 6.190253e-08
#> 802       TG 20:4_33:2          Triacylglycerols -1.3798243260 4.999836e-09
#> 803       TG 20:4_34:0          Triacylglycerols -0.8973961559 1.942399e-08
#> 804       TG 20:4_34:1          Triacylglycerols -1.5574121489 1.279990e-10
#> 805       TG 20:4_34:2          Triacylglycerols -1.7696217886 6.291336e-10
#> 806       TG 20:4_34:3          Triacylglycerols -1.7062178661 1.587741e-09
#> 807       TG 20:4_35:3          Triacylglycerols -1.0281296368 4.684362e-06
#> 808       TG 20:4_36:4          Triacylglycerols -0.9924061052 1.917087e-07
#> 809       TG 20:4_36:5          Triacylglycerols -1.1845389197 4.501141e-09
#> 810       TG 20:5_34:1          Triacylglycerols -1.5094443227 1.736480e-09
#> 811       TG 20:5_36:2          Triacylglycerols -1.6602411269 1.337116e-09
#> 812       TG 20:5_36:3          Triacylglycerols -1.6347576000 2.308948e-08
#> 813       TG 22:0_32:4          Triacylglycerols -0.8557222679 1.343958e-05
#> 814       TG 22:1_32:5          Triacylglycerols -1.3147235625 2.447208e-06
#> 815       TG 22:2_32:4          Triacylglycerols -0.9588565744 2.822861e-03
#> 816       TG 22:3_30:2          Triacylglycerols -1.0980133938 5.707067e-05
#> 817       TG 22:4_32:0          Triacylglycerols -1.1787544334 2.746774e-06
#> 818       TG 22:4_32:2          Triacylglycerols -1.1157332999 7.845193e-09
#> 819       TG 22:4_34:2          Triacylglycerols -1.5627314431 4.208208e-09
#> 820       TG 22:5_32:0          Triacylglycerols -0.8113686162 7.266361e-08
#> 821       TG 22:5_32:1          Triacylglycerols -1.3534494207 2.971020e-10
#> 822       TG 22:5_34:1          Triacylglycerols -1.5928490934 2.907066e-15
#> 823       TG 22:6_32:0          Triacylglycerols -1.1194768245 1.557580e-08
#> 824       TG 22:6_32:1          Triacylglycerols -1.4536445968 1.032540e-08
#> 825       TG 22:6_34:1          Triacylglycerols -1.5552854227 8.849741e-12
#> 826       TG 22:6_34:2          Triacylglycerols -1.6848124159 9.100375e-10
#> 827                Thr                Aminoacids  1.2180279083 5.522918e-15
#> 828               TMAO              Amine Oxides -1.0349951867 5.800442e-03
#> 829               TMCA                Bile Acids -1.2669195327 2.236198e-04
#> 830       Trigonelline                 Alkaloids  0.0758405622 7.358387e-01
#> 831                Trp                Aminoacids -0.8414685883 1.132489e-11
#> 832         TrpBetaine        Aminoacids Related -1.0362066129 7.745687e-05
#> 833                Tyr                Aminoacids  0.1589710375 1.006383e-02
#> 834                Val                Aminoacids  0.1202717820 7.484468e-02
#>             qval
#> 1   1.344839e-04
#> 2   8.821852e-10
#> 3   1.033145e-05
#> 4   1.281497e-05
#> 5   1.544771e-15
#> 6   1.367397e-09
#> 7   2.916968e-01
#> 8   8.310698e-09
#> 9   9.394881e-01
#> 10  1.218615e-09
#> 11  1.938536e-11
#> 12  1.705170e-03
#> 13  6.882992e-02
#> 14  9.935905e-04
#> 15  1.455796e-15
#> 16  3.411643e-15
#> 17  3.042986e-02
#> 18  1.188642e-08
#> 19  2.763212e-04
#> 20  2.366136e-03
#> 21  8.367249e-04
#> 22  6.934566e-05
#> 23  6.324912e-01
#> 24  2.172568e-04
#> 25  1.185260e-02
#> 26  6.535440e-05
#> 27  2.993578e-01
#> 28  4.528047e-02
#> 29  4.655610e-02
#> 30  7.147421e-02
#> 31  2.020926e-01
#> 32  1.078022e-02
#> 33  9.198232e-04
#> 34  6.562698e-01
#> 35  4.904602e-01
#> 36  1.415972e-02
#> 37  4.273825e-12
#> 38  2.724751e-03
#> 39  2.997158e-04
#> 40  3.601756e-07
#> 41  4.590153e-12
#> 42  7.409791e-01
#> 43  3.750989e-12
#> 44  7.490850e-01
#> 45  5.383169e-02
#> 46  2.812025e-03
#> 47  9.485924e-01
#> 48  6.417034e-04
#> 49  2.766109e-07
#> 50  2.987733e-02
#> 51  4.660563e-02
#> 52  6.606251e-02
#> 53  6.618333e-02
#> 54  1.113978e-09
#> 55  1.145417e-03
#> 56  8.938667e-06
#> 57  1.217628e-01
#> 58  3.068137e-03
#> 59  2.831123e-04
#> 60  2.196770e-06
#> 61  8.662034e-04
#> 62  8.912076e-04
#> 63  3.596710e-04
#> 64  3.821502e-02
#> 65  6.463683e-07
#> 66  1.251395e-12
#> 67  2.513442e-14
#> 68  2.419512e-03
#> 69  2.688598e-09
#> 70  7.088317e-12
#> 71  4.451780e-08
#> 72  5.787897e-04
#> 73  4.132699e-02
#> 74  7.175416e-01
#> 75  8.665783e-02
#> 76  1.128783e-04
#> 77  9.464375e-04
#> 78  1.139497e-06
#> 79  1.201582e-01
#> 80  6.361099e-04
#> 81  2.917460e-08
#> 82  6.366410e-01
#> 83  9.075413e-07
#> 84  2.650505e-07
#> 85  3.931010e-04
#> 86  4.523508e-08
#> 87  5.451466e-01
#> 88  1.611972e-12
#> 89  7.947284e-15
#> 90  1.110425e-13
#> 91  1.817055e-02
#> 92  8.376470e-09
#> 93  2.291252e-14
#> 94  4.203153e-12
#> 95  3.811520e-14
#> 96  3.805543e-11
#> 97  1.574291e-14
#> 98  7.387964e-15
#> 99  4.651926e-12
#> 100 3.953442e-15
#> 101 7.460242e-15
#> 102 5.800082e-04
#> 103 3.811520e-14
#> 104 2.706273e-13
#> 105 8.018795e-01
#> 106 4.050899e-10
#> 107 6.176314e-13
#> 108 8.514834e-07
#> 109 3.600752e-14
#> 110 2.723213e-13
#> 111 6.336782e-08
#> 112 4.949609e-10
#> 113 2.408318e-01
#> 114 6.965740e-01
#> 115 9.032755e-09
#> 116 1.217628e-01
#> 117 9.108681e-03
#> 118 1.759786e-05
#> 119 6.103640e-04
#> 120 8.923101e-02
#> 121 1.328175e-08
#> 122 3.732717e-06
#> 123 1.258340e-01
#> 124 4.973746e-09
#> 125 8.466454e-01
#> 126 4.612777e-11
#> 127 1.510233e-11
#> 128 2.653635e-11
#> 129 2.206048e-01
#> 130 8.720812e-01
#> 131 3.737001e-08
#> 132 1.753997e-08
#> 133 7.106083e-03
#> 134 2.163337e-02
#> 135 2.225722e-15
#> 136 7.763888e-01
#> 137 8.807833e-14
#> 138 5.490500e-01
#> 139 8.154690e-10
#> 140 4.740395e-06
#> 141 2.389821e-08
#> 142 4.351819e-01
#> 143 1.063580e-02
#> 144 9.261550e-09
#> 145 1.157437e-06
#> 146 1.062738e-01
#> 147 9.094225e-01
#> 148 3.872000e-03
#> 149 1.211938e-01
#> 150 1.988450e-01
#> 151 6.385303e-02
#> 152 7.942913e-01
#> 153 9.505075e-01
#> 154 1.259311e-01
#> 155 4.324598e-03
#> 156 1.124170e-07
#> 157 2.088485e-01
#> 158 5.767030e-01
#> 159 1.129048e-01
#> 160 2.747070e-01
#> 161 3.545401e-01
#> 162 9.102488e-01
#> 163 1.964395e-07
#> 164 1.653214e-07
#> 165 1.490458e-02
#> 166 2.034315e-02
#> 167 7.553736e-01
#> 168 2.472099e-06
#> 169 1.841823e-06
#> 170 1.849297e-06
#> 171 1.784275e-15
#> 172 3.259033e-07
#> 173 2.677595e-14
#> 174 1.447608e-06
#> 175 3.517445e-04
#> 176 6.776958e-06
#> 177 6.683708e-09
#> 178 4.611007e-16
#> 179 1.259032e-13
#> 180 4.290626e-05
#> 181 3.108825e-08
#> 182 1.655135e-11
#> 183 1.839410e-02
#> 184 5.560833e-06
#> 185 1.126089e-13
#> 186 4.593046e-15
#> 187 4.777923e-13
#> 188 5.909480e-05
#> 189 6.734178e-14
#> 190 3.854845e-14
#> 191 9.772164e-13
#> 192 5.480528e-15
#> 193 1.343860e-14
#> 194 3.017413e-15
#> 195 3.244042e-15
#> 196 1.519160e-14
#> 197 2.758032e-07
#> 198 1.147649e-07
#> 199 8.356817e-11
#> 200 1.326309e-04
#> 201 2.349994e-13
#> 202 5.787897e-04
#> 203 8.260589e-13
#> 204 1.615025e-12
#> 205 3.634231e-14
#> 206 5.572998e-14
#> 207 3.354578e-15
#> 208 1.382108e-14
#> 209 8.979065e-15
#> 210 1.823917e-14
#> 211 5.324511e-07
#> 212 1.136927e-04
#> 213 3.332202e-01
#> 214 4.392717e-06
#> 215 3.273130e-06
#> 216 1.128081e-05
#> 217 2.789118e-04
#> 218 1.912245e-05
#> 219 7.044562e-03
#> 220 7.366454e-05
#> 221 4.163346e-02
#> 222 2.242355e-09
#> 223 9.937353e-07
#> 224 8.309764e-01
#> 225 2.815471e-13
#> 226 4.193770e-08
#> 227 2.866840e-04
#> 228 4.833756e-06
#> 229 1.479578e-13
#> 230 3.126037e-01
#> 231 1.443759e-06
#> 232 1.452986e-12
#> 233 1.000918e-07
#> 234 3.116244e-12
#> 235 1.552603e-11
#> 236 1.004786e-12
#> 237 2.891088e-05
#> 238 2.554371e-06
#> 239 5.124745e-03
#> 240 8.832881e-02
#> 241 4.679784e-08
#> 242 8.919955e-09
#> 243 5.806845e-12
#> 244 1.784275e-15
#> 245 3.477058e-13
#> 246 8.151933e-11
#> 247 2.242585e-13
#> 248 1.038198e-13
#> 249 3.332393e-09
#> 250 1.099352e-01
#> 251 3.292170e-08
#> 252 4.538198e-09
#> 253 2.545766e-08
#> 254 3.350079e-11
#> 255 7.119398e-10
#> 256 1.584868e-14
#> 257 1.365554e-06
#> 258 1.759786e-05
#> 259 2.889361e-07
#> 260 5.460654e-10
#> 261 6.440505e-12
#> 262 3.090684e-13
#> 263 1.735459e-06
#> 264 8.099953e-11
#> 265 1.544771e-15
#> 266 3.256864e-15
#> 267 1.583230e-15
#> 268 1.627564e-14
#> 269 5.667550e-07
#> 270 3.860605e-09
#> 271 1.455796e-15
#> 272 2.138241e-14
#> 273 2.520863e-07
#> 274 5.027917e-15
#> 275 7.066133e-16
#> 276 2.168191e-15
#> 277 5.592685e-02
#> 278 4.196946e-02
#> 279 8.570652e-12
#> 280 1.369124e-10
#> 281 1.544771e-15
#> 282 2.647345e-07
#> 283 4.589642e-05
#> 284 5.976893e-04
#> 285 1.305558e-11
#> 286 4.549811e-07
#> 287 8.118544e-08
#> 288 5.626250e-04
#> 289 2.186373e-01
#> 290 2.528752e-01
#> 291 1.372837e-08
#> 292 1.287258e-09
#> 293 8.533638e-01
#> 294 1.515900e-01
#> 295 3.244042e-15
#> 296 5.787303e-12
#> 297 1.891905e-06
#> 298 1.826359e-12
#> 299 6.248355e-03
#> 300 2.168191e-15
#> 301 1.801268e-13
#> 302 5.051701e-02
#> 303 9.431149e-15
#> 304 3.927253e-15
#> 305 1.801424e-09
#> 306 3.515115e-02
#> 307 2.509513e-14
#> 308 3.092172e-06
#> 309 4.271723e-04
#> 310 2.186373e-01
#> 311 5.108460e-11
#> 312 9.251660e-02
#> 313 3.984853e-03
#> 314 1.918954e-07
#> 315 1.173255e-11
#> 316 3.133927e-01
#> 317 1.131764e-02
#> 318 7.820866e-02
#> 319 2.257753e-15
#> 320 9.440462e-10
#> 321 1.627564e-14
#> 322 2.168191e-15
#> 323 5.015589e-15
#> 324 1.739042e-12
#> 325 5.169123e-09
#> 326 1.376618e-03
#> 327 1.188799e-03
#> 328 4.110427e-06
#> 329 5.332081e-07
#> 330 8.631372e-10
#> 331 7.273220e-01
#> 332 2.280530e-09
#> 333 3.862661e-13
#> 334 1.600449e-07
#> 335 3.279558e-05
#> 336 1.854978e-03
#> 337 1.556944e-09
#> 338 2.606075e-07
#> 339 1.027176e-10
#> 340 2.131138e-10
#> 341 1.105459e-08
#> 342 3.696567e-13
#> 343 2.257753e-15
#> 344 2.067074e-11
#> 345 3.510063e-02
#> 346 8.794819e-04
#> 347 1.523564e-05
#> 348 3.866707e-02
#> 349 1.574291e-14
#> 350 2.089201e-14
#> 351 4.343448e-09
#> 352 7.685387e-05
#> 353 1.561802e-05
#> 354 1.077003e-14
#> 355 1.540662e-08
#> 356 1.544771e-15
#> 357 4.597870e-15
#> 358 2.056753e-11
#> 359 4.119144e-15
#> 360 1.023992e-14
#> 361 1.810445e-14
#> 362 1.003611e-14
#> 363 9.715568e-12
#> 364 5.044414e-02
#> 365 3.085767e-11
#> 366 1.301684e-13
#> 367 7.040989e-10
#> 368 2.066851e-05
#> 369 1.760470e-07
#> 370 4.578691e-06
#> 371 6.739142e-05
#> 372 6.446648e-15
#> 373 1.040054e-03
#> 374 4.713935e-09
#> 375 2.121437e-02
#> 376 2.497025e-01
#> 377 6.723922e-03
#> 378 5.295737e-04
#> 379 4.282517e-05
#> 380 1.077092e-01
#> 381 1.396687e-01
#> 382 2.205921e-02
#> 383 4.720518e-06
#> 384 1.160403e-08
#> 385 4.456773e-03
#> 386 8.875991e-02
#> 387 1.726310e-01
#> 388 8.832207e-01
#> 389 1.671083e-06
#> 390 2.483546e-03
#> 391 2.090233e-11
#> 392 1.142425e-05
#> 393 9.421009e-05
#> 394 1.522352e-15
#> 395 1.801424e-09
#> 396 1.583230e-15
#> 397 2.029923e-15
#> 398 1.930416e-15
#> 399 2.168191e-15
#> 400 1.930416e-15
#> 401 2.354093e-15
#> 402 2.168191e-15
#> 403 1.291942e-15
#> 404 1.544771e-15
#> 405 3.696567e-13
#> 406 1.182219e-15
#> 407 1.182219e-15
#> 408 1.182219e-15
#> 409 1.455796e-15
#> 410 1.182219e-15
#> 411 1.182219e-15
#> 412 3.398323e-15
#> 413 4.062294e-15
#> 414 6.420480e-16
#> 415 5.948271e-16
#> 416 7.066133e-16
#> 417 1.455796e-15
#> 418 1.182219e-15
#> 419 2.029923e-15
#> 420 6.711198e-14
#> 421 3.244042e-15
#> 422 1.376631e-15
#> 423 1.115044e-15
#> 424 1.182219e-15
#> 425 1.182219e-15
#> 426 1.480736e-14
#> 427 1.241759e-15
#> 428 4.656748e-15
#> 429 2.029923e-15
#> 430 2.248136e-15
#> 431 1.455796e-15
#> 432 1.930416e-15
#> 433 1.546811e-15
#> 434 1.241759e-15
#> 435 2.729209e-15
#> 436 1.241759e-15
#> 437 4.593046e-15
#> 438 2.129761e-14
#> 439 5.600036e-14
#> 440 2.488598e-12
#> 441 2.985219e-15
#> 442 4.593046e-15
#> 443 1.455796e-15
#> 444 1.544933e-14
#> 445 2.521513e-15
#> 446 1.100359e-13
#> 447 9.531319e-15
#> 448 2.100606e-15
#> 449 2.100606e-15
#> 450 1.544771e-15
#> 451 1.875326e-15
#> 452 1.232705e-15
#> 453 3.498609e-15
#> 454 1.455796e-15
#> 455 4.465537e-15
#> 456 9.707240e-14
#> 457 1.343860e-14
#> 458 2.310916e-04
#> 459 2.665445e-14
#> 460 1.008105e-04
#> 461 7.477256e-14
#> 462 1.130533e-14
#> 463 1.462242e-14
#> 464 1.339229e-14
#> 465 7.248804e-12
#> 466 5.153428e-10
#> 467 8.193501e-16
#> 468 1.115044e-15
#> 469 1.182219e-15
#> 470 2.211418e-12
#> 471 2.248136e-15
#> 472 2.277728e-12
#> 473 5.979329e-14
#> 474 2.057937e-14
#> 475 3.996206e-14
#> 476 1.938232e-01
#> 477 2.355369e-10
#> 478 1.706837e-12
#> 479 4.510641e-05
#> 480 1.547005e-11
#> 481 1.097312e-04
#> 482 1.524783e-01
#> 483 3.935796e-01
#> 484 9.410176e-09
#> 485 1.017954e-14
#> 486 1.381978e-11
#> 487 4.119144e-15
#> 488 2.729209e-15
#> 489 7.129204e-14
#> 490 9.902368e-01
#> 491 1.455796e-15
#> 492 9.622362e-15
#> 493 1.701249e-14
#> 494 1.403540e-14
#> 495 2.057799e-03
#> 496 1.084532e-13
#> 497 9.024186e-12
#> 498 4.593046e-15
#> 499 1.773952e-15
#> 500 2.614568e-04
#> 501 1.801298e-13
#> 502 6.420480e-16
#> 503 4.465537e-15
#> 504 1.225315e-12
#> 505 1.527433e-15
#> 506 3.810479e-14
#> 507 5.071308e-11
#> 508 7.457855e-08
#> 509 6.791744e-03
#> 510 2.029923e-15
#> 511 2.418119e-10
#> 512 4.117521e-13
#> 513 1.870365e-13
#> 514 2.770669e-11
#> 515 8.425559e-10
#> 516 6.971494e-08
#> 517 1.286373e-14
#> 518 4.103624e-14
#> 519 1.068497e-14
#> 520 7.066133e-16
#> 521 5.090428e-14
#> 522 6.868742e-11
#> 523 5.882458e-13
#> 524 1.494413e-06
#> 525 7.223266e-06
#> 526 9.372049e-02
#> 527 3.649023e-06
#> 528 7.276702e-13
#> 529 2.327115e-13
#> 530 3.094695e-01
#> 531 2.514345e-04
#> 532 1.395629e-03
#> 533 1.620122e-11
#> 534 3.477058e-13
#> 535 2.045808e-10
#> 536 4.156172e-13
#> 537 9.236386e-02
#> 538 4.029527e-04
#> 539 1.526386e-11
#> 540 1.669845e-05
#> 541 5.900296e-09
#> 542 2.522269e-06
#> 543 1.825618e-05
#> 544 2.383965e-01
#> 545 4.597870e-15
#> 546 2.765202e-02
#> 547 1.471963e-06
#> 548 5.755304e-11
#> 549 8.850819e-14
#> 550 2.565040e-09
#> 551 1.875326e-15
#> 552 1.996896e-07
#> 553 1.455796e-15
#> 554 1.724122e-10
#> 555 5.468350e-15
#> 556 4.840138e-05
#> 557 1.374879e-14
#> 558 1.401896e-12
#> 559 6.170734e-08
#> 560 1.814375e-14
#> 561 2.729209e-15
#> 562 1.638632e-04
#> 563 2.615559e-10
#> 564 1.714364e-13
#> 565 1.203215e-10
#> 566 1.431383e-14
#> 567 2.667679e-07
#> 568 7.894246e-14
#> 569 9.496394e-02
#> 570 7.092678e-12
#> 571 4.333615e-08
#> 572 6.858928e-10
#> 573 2.427006e-15
#> 574 6.037226e-15
#> 575 1.317788e-11
#> 576 1.509281e-11
#> 577 1.144350e-14
#> 578 4.878092e-15
#> 579 3.838540e-12
#> 580 4.119144e-15
#> 581 2.211426e-11
#> 582 4.119144e-15
#> 583 1.826359e-12
#> 584 5.606622e-14
#> 585 1.752311e-09
#> 586 3.322097e-11
#> 587 1.455796e-15
#> 588 1.624156e-01
#> 589 1.271240e-03
#> 590 3.214774e-09
#> 591 1.656232e-07
#> 592 1.463318e-04
#> 593 5.017559e-07
#> 594 8.696421e-04
#> 595 7.851263e-03
#> 596 9.890540e-14
#> 597 1.424261e-15
#> 598 1.544771e-15
#> 599 1.182219e-15
#> 600 1.232705e-15
#> 601 1.455796e-15
#> 602 5.848627e-15
#> 603 5.009986e-04
#> 604 7.066133e-16
#> 605 1.130533e-14
#> 606 1.013263e-14
#> 607 2.363327e-07
#> 608 3.138032e-14
#> 609 9.003658e-15
#> 610 1.055693e-14
#> 611 7.153630e-14
#> 612 1.691107e-11
#> 613 4.569422e-05
#> 614 5.226148e-08
#> 615 2.502389e-12
#> 616 6.051965e-11
#> 617 7.181121e-08
#> 618 2.284928e-04
#> 619 1.592775e-09
#> 620 3.670265e-03
#> 621 4.968704e-06
#> 622 3.339401e-01
#> 623 3.095232e-06
#> 624 6.143434e-02
#> 625 7.378636e-01
#> 626 2.728613e-09
#> 627 3.931010e-04
#> 628 1.662554e-02
#> 629 2.677417e-09
#> 630 1.097420e-11
#> 631 4.196104e-09
#> 632 5.949906e-01
#> 633 3.296880e-01
#> 634 1.154993e-06
#> 635 7.086867e-15
#> 636 1.100359e-13
#> 637 2.060267e-06
#> 638 2.109829e-15
#> 639 9.727920e-04
#> 640 3.408802e-04
#> 641 3.149088e-07
#> 642 2.737925e-04
#> 643 9.418622e-07
#> 644 2.183166e-09
#> 645 3.920097e-09
#> 646 5.963870e-03
#> 647 4.062558e-06
#> 648 1.094139e-04
#> 649 1.230079e-11
#> 650 3.156471e-11
#> 651 7.235888e-09
#> 652 6.882271e-08
#> 653 1.679813e-01
#> 654 6.252695e-05
#> 655 5.944625e-05
#> 656 3.155542e-07
#> 657 2.024417e-01
#> 658 8.353579e-05
#> 659 6.120456e-08
#> 660 3.243213e-08
#> 661 1.032563e-01
#> 662 1.302467e-05
#> 663 8.358256e-03
#> 664 8.922171e-07
#> 665 1.017328e-09
#> 666 1.066709e-09
#> 667 1.992082e-08
#> 668 1.560327e-05
#> 669 6.023412e-10
#> 670 8.444902e-09
#> 671 1.087373e-10
#> 672 3.731809e-10
#> 673 2.109418e-07
#> 674 1.634536e-06
#> 675 5.332150e-05
#> 676 8.201406e-10
#> 677 1.173255e-11
#> 678 1.384803e-11
#> 679 8.360144e-11
#> 680 7.582215e-11
#> 681 1.083947e-08
#> 682 1.985801e-11
#> 683 7.713060e-14
#> 684 3.135811e-03
#> 685 2.702045e-04
#> 686 2.712308e-01
#> 687 1.303586e-04
#> 688 3.197725e-08
#> 689 1.565142e-01
#> 690 5.665685e-02
#> 691 1.724924e-08
#> 692 8.100719e-10
#> 693 3.825229e-09
#> 694 1.645583e-06
#> 695 8.653144e-12
#> 696 2.502304e-09
#> 697 2.418119e-10
#> 698 7.937541e-10
#> 699 2.669482e-02
#> 700 1.063137e-08
#> 701 1.791453e-04
#> 702 6.015404e-09
#> 703 1.232696e-11
#> 704 4.053761e-09
#> 705 3.027755e-12
#> 706 1.415285e-05
#> 707 4.914196e-04
#> 708 5.046616e-03
#> 709 4.755077e-06
#> 710 1.935707e-06
#> 711 4.468785e-05
#> 712 1.322209e-07
#> 713 1.377274e-05
#> 714 6.432690e-02
#> 715 3.113742e-01
#> 716 1.893010e-02
#> 717 3.369262e-04
#> 718 1.370803e-02
#> 719 5.497477e-04
#> 720 3.790285e-06
#> 721 5.046082e-08
#> 722 8.679521e-09
#> 723 2.192863e-03
#> 724 3.627118e-07
#> 725 4.571946e-11
#> 726 8.821852e-10
#> 727 1.328175e-08
#> 728 1.712972e-10
#> 729 2.429572e-05
#> 730 7.284808e-07
#> 731 5.922048e-07
#> 732 1.394442e-09
#> 733 1.196656e-09
#> 734 8.968083e-04
#> 735 3.172323e-07
#> 736 3.307979e-09
#> 737 2.719431e-12
#> 738 9.970250e-11
#> 739 7.105918e-07
#> 740 1.899582e-10
#> 741 1.927211e-11
#> 742 8.117715e-09
#> 743 3.027755e-12
#> 744 6.711198e-14
#> 745 3.060829e-13
#> 746 5.446104e-11
#> 747 3.547554e-10
#> 748 3.197917e-05
#> 749 6.019199e-07
#> 750 1.461797e-11
#> 751 2.874992e-08
#> 752 1.023992e-14
#> 753 7.153630e-14
#> 754 1.629352e-07
#> 755 2.211256e-07
#> 756 9.068614e-12
#> 757 2.552475e-10
#> 758 3.705146e-07
#> 759 2.184301e-09
#> 760 1.112447e-08
#> 761 2.520863e-07
#> 762 2.116331e-12
#> 763 1.471869e-09
#> 764 2.101619e-09
#> 765 2.739856e-09
#> 766 3.212268e-08
#> 767 6.703983e-07
#> 768 8.722855e-10
#> 769 3.207174e-10
#> 770 1.283869e-10
#> 771 2.696731e-07
#> 772 9.715271e-08
#> 773 1.017328e-09
#> 774 5.459555e-08
#> 775 9.589520e-09
#> 776 5.382603e-06
#> 777 2.269624e-04
#> 778 2.308017e-05
#> 779 8.912076e-04
#> 780 4.833756e-06
#> 781 5.154846e-07
#> 782 1.704823e-09
#> 783 6.419473e-10
#> 784 8.784318e-11
#> 785 2.904691e-05
#> 786 4.938803e-10
#> 787 3.862661e-13
#> 788 1.903591e-10
#> 789 1.152760e-05
#> 790 7.959979e-08
#> 791 3.632090e-10
#> 792 3.836236e-09
#> 793 3.578505e-01
#> 794 2.225786e-10
#> 795 2.942552e-11
#> 796 9.363172e-12
#> 797 7.109688e-10
#> 798 1.179983e-08
#> 799 2.524474e-06
#> 800 2.501058e-06
#> 801 1.086878e-07
#> 802 9.788411e-09
#> 803 3.576073e-08
#> 804 3.167691e-10
#> 805 1.406695e-09
#> 806 3.335457e-09
#> 807 6.853961e-06
#> 808 3.172323e-07
#> 809 8.919955e-09
#> 810 3.638754e-09
#> 811 2.837545e-09
#> 812 4.222945e-08
#> 813 1.906226e-05
#> 814 3.690726e-06
#> 815 3.411979e-03
#> 816 7.764591e-05
#> 817 4.110427e-06
#> 818 1.497229e-08
#> 819 8.416415e-09
#> 820 1.267813e-07
#> 821 7.059347e-10
#> 822 1.584868e-14
#> 823 2.912604e-08
#> 824 1.952694e-08
#> 825 2.501927e-11
#> 826 1.992051e-09
#> 827 2.725511e-14
#> 828 6.891123e-03
#> 829 2.914045e-04
#> 830 7.493156e-01
#> 831 3.148319e-11
#> 832 1.046986e-04
#> 833 1.180483e-02
#> 834 8.311647e-02
```
