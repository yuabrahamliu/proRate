FILE=config
RSCRIPT=`cat $(FILE) | grep 'Rscript' | cut -d ' ' -f 3`

.PONEY: all

all:
	$(RSCRIPT) elongHMM.R