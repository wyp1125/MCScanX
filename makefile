CC = g++
CFLAGS =

prefix = /usr/local
bindir = $(prefix)/bin

all: MCScanX duplicate_gene_classifier downstream_analyses
	if [ ! -d ./bin ]; then mkdir -p ./bin; fi
	cp -f MCScanX duplicate_gene_classifier ./bin

MCScanX: struct.o mcscan.o read_data.o out_utils.o dagchainer.o msa.o permutation.o
	$(CC) $(CFLAGS) struct.o mcscan.o read_data.o out_utils.o dagchainer.o msa.o permutation.o -o MCScanX

duplicate_gene_classifier: struct.o dup_classifier.o read_data.o out_utils.o dagchainer.o cls.o permutation.o
	$(CC) $(CFLAGS) struct.o dup_classifier.o read_data.o out_utils.o dagchainer.o cls.o permutation.o -o duplicate_gene_classifier

downstream_analyses: dissect_multiple_alignment.o detect_syntenic_tandem_arrays.o
	cd downstream_analyses/ && ${MAKE}

clean:
	rm -f *.o MCScanX duplicate_gene_classifier
	cd downstream_analyses/ && ${MAKE} clean

install: MCScanX duplicate_gene_classifier downstream_analyses
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	cp -f MCScanX duplicate_gene_classifier $(bindir)
	cd downstream_analyses/ && ${MAKE} install

uninstall:
	cd $(bindir); \
	rm -f MCScanX duplicate_gene_classifier
	cd downstream_analyses/ && ${MAKE} uninstall