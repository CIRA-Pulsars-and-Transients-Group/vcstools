TARGETS = vcstools-tree.pdf vcstools-tree.png

all: $(TARGETS)

vcstools-tree.pdf: vcstools-tree.tex
	pdflatex $<
	pdflatex $<

vcstools-tree.png: vcstools-tree.pdf
	convert -flatten -background white -density 400 $< $@

clean:
	$(RM) $(TARGETS)
	$(RM) *.aux *.log *synctex*
