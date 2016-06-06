PDFS=sources.pdf

all: $(PDFS)

%.pdf: %.md
	pandoc -o $@ $<

.PHONY: clean
clean:
	rm -f $(PDFS)
