.PHONY: all
all: k2l example

.PHONY: k2l
k2l: 
	(cd ./src; $(MAKE))

.PHONY: example
example: 
	(cd ./example; $(MAKE))

.PHONY: clean
clean: cleank2l cleanexample

.PHONY: cleank2l
cleank2l:
	(cd ./src; $(MAKE) clean)

.PHONY: cleanexample
cleanexample:
	(cd ./example; $(MAKE) clean)		
