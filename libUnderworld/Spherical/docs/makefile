INPUTFILE=rotatingBC.markdown
OUTPUTPDF=doc.pdf
OUTPUTHTML=doc.html

OPTION=
HTMLOPTION:=$(OPTION) --webtex=http://latex.codecogs.com/gif.latex?


all: $(INPUTFILE)
	pandoc $(INPUTFILE) -o $(OUTPUTPDF) $(OPTION)
	pandoc $(INPUTFILE) -o $(OUTPUTHTML) $(HTMLOPTION)

clean:
	rm $(OUTPUTFILE)
