all: freesurfer_r_journal_submission.zip
freesurfer_r_journal_submission.zip: \
	RJwrapper.pdf \
	RJwrapper.tex \
	muschelli.Rmd \
	muschelli.bib 
	R -e "rmarkdown::render('muschelli.Rmd')"	
	R -e "knitr::purl('muschelli.Rmd')"
	zip -r freesurfer_r_journal_submission.zip \
	RJwrapper.pdf RJwrapper.tex muschelli.tex \
	mask.nii.gz muschelli.bib muschelli.Rmd \
	muschelli.R muschelli_files/*


clean: 
	rm freesurfer_r_journal_submission.zip
