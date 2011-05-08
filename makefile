mcscan2: struct.cc mcscan.cc read_data.cc out_utils.cc dagchainer.cc msa.cc permutation.cc	
	g++ struct.cc mcscan.cc read_data.cc out_utils.cc dagchainer.cc msa.cc permutation.cc -o MCScanX
	g++ struct.cc dup_classifier.cc read_data.cc out_utils.cc dagchainer.cc cls.cc permutation.cc -o duplicate_gene_classifier
	g++ dissect_multiple_alignment.cc -o downstream_analyses/dissect_multiple_alignment
	g++ detect_syntenic_tandem_arrays.cc -o downstream_analyses/detect_syntenic_tandem_arrays
	cd downstream_analyses/ && ${MAKE}

