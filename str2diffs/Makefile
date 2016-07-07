

RMDFILE=README


all:
	@echo "Spinning..."
	Rscript -e "knitr::spin('$(RMDFILE).R', knit = FALSE)"
	Rscript -e "rmarkdown::render('$(RMDFILE).Rmd')"

clean:
	@echo "Cleaning..."
	rm -f *.Rmd *.html *.md
