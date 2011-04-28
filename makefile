mcscan2: struct.cc mcscan.cc read_data.cc out_utils.cc dagchainer.cc msa.cc permutation.cc	
	g++ struct.cc mcscan.cc read_data.cc out_utils.cc dagchainer.cc msa.cc permutation.cc -o bin/mcscan2
	g++ struct.cc dup_classifier.cc read_data.cc out_utils.cc dagchainer.cc cls.cc permutation.cc -o bin/duplicate_classifier
	g++ species_specific_synteny.cc -o post_processing/species_specific_synteny
	g++ species_specific_tandem.cc -o post_processing/species_specific_tandem
