all: 
	(cd src; $(MAKE))
	(cd example; $(MAKE))

.PHONY: all clean
clean:
	(cd src; $(MAKE) clean)
	(cd example; $(MAKE) clean)
	