CC = g++
CFLAGS =
JFLAGS = -g
JC = javac
.SUFFIXES: .java .class

prefix = /usr/local
bindir = $(prefix)/bin

all: classes dissect_multiple_alignment detect_syntenic_tandem_arrays
	if [ ! -d ./bin ]; then mkdir -p ../bin; fi
	cp -f dissect_multiple_alignment detect_syntenic_tandem_arrays *.class *.pl ../bin

dissect_multiple_alignment: ../dissect_multiple_alignment.o
	$(CC) $(CFLAGS) ../dissect_multiple_alignment.o -o dissect_multiple_alignment

detect_syntenic_tandem_arrays: ../detect_syntenic_tandem_arrays.o
	$(CC) $(CFLAGS) ../detect_syntenic_tandem_arrays.o -o detect_syntenic_tandem_arrays

CLASSES = $(wildcard *.java)
classes: $(CLASSES:.java=.class)
.java.class:
	$(JC) $(JFLAGS) $*.java

clean:
	$(RM) *.class *.o dissect_multiple_alignment detect_syntenic_tandem_arrays

install: dissect_multiple_alignment detect_syntenic_tandem_arrays classes *.pl
	cp -f dissect_multiple_alignment detect_syntenic_tandem_arrays *.class *.pl $(bindir)

uninstall:
	cd $(bindir); \
	rm -f dissect_multiple_alignment detect_syntenic_tandem_arrays add_kaks_to_synteny.pl detect_synteny_within_gene_families.pl group_syntenic_genes.pl origin_enrichment_analysis.pl bar_plotter.class circle_plotter.class Cubic.class dot_plotter.class dual_synteny_plotter.class family_circle_plotter.class family_tree_plotter.class