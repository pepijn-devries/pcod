# make the outputs from the elicitation into a proper .rda object
rm(list = ls())
load("../data-raw/HS_mature_Days.robj")
hsmaturedays <- temp
save(hsmaturedays, file = '../data/hsMaturedays.rda')

rm(list = ls())
load("../data-raw/HSjuve_Days.robj") 
hsjuvedays <- temp
save(hsjuvedays, file = '../data/hsJuvedays.rda')

rm(list = ls())
load("../data-raw/HSdependent_Days.robj") 
hsdependentdays <- temp
save(hsdependentdays, file = '../data/hsDependentdays.rda')


rm(list = ls())
load("../data-raw/output_matureHS_290113.robj")
hsoutput_mature <- output_HarbourSeal_mature
save(hsoutput_mature, file = '../data/output_HarbourSeal_mature.rda')

rm(list = ls())
load("../data-raw/output_juvenileHS_290113.robj")                   
hsoutput_juvenile <- output_HarbourSeal_juvenile
save(hsoutput_juvenile, file = '../data/output_HarbourSeal_juvenile.rda')

rm(list = ls())
load("../data-raw/output_dependentHS_290113.robj") 
hsoutput_dependent <- output_HarbourSeal_dependent
save(hsoutput_dependent, file = '../data/output_HarbourSeal_dependent.rda')
