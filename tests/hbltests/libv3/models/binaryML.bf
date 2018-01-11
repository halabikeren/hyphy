LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/stats.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/mpi.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("libv3/models/binary/charbinary.bf");
LoadFunctionLibrary("libv3/models/binary.bf");
LoadFunctionLibrary("libv3/models/binary/empirical.bf");

function test_binaryML() {

	// Test dataset
	data_file_name = "`PATH_TO_CURRENT_BF`../data/binary.fas";
	tree_file_name = "`PATH_TO_CURRENT_BF`../data/mammals.tre";

    // Read in the nucleotide alignments into hky85.nuc_data and hk85.nuc_filter variables
    binary_character_info = alignments.ReadNucleotideAlignment(data_file_name, "binaryML.char_data", "binaryML.char_filter");
    char_data_default = { "0" : "binaryML.char_filter"};

    name_space = & model_binaryML;

    // define the model
    binaryML.model = model.generic.DefineModel("models.binaryML.ModelDescription",
        name_space, {
            "0": terms.local
        },
        char_data_default,
        None);

    //fprintf(stdout, hky85.model);
	utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
	ExecuteCommands ('binaryML.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions(binary_character_info[terms.data.partitions], binaryML.name_mapping);',
						 {"0" : tree_file_name});                    
	utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);
    binaryML.trees = utility.Map(binaryML.partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');

    lf_components = {1,2};
    lf_components[0] = "binaryML.char_filter";
    lf_components[1] = "tree_0";

    tree = binaryML.trees["0"];

    model_assignment = {
        terms.default: binaryML.model
    };

    model.ApplyModelToTree(lf_components[1], 
                           tree, 
                           model_assignment, 
                           None);

    model_id_to_object = {
        name_space: binaryML.model
    };

    LikelihoodFunction likelihoodFunction = (lf_components);
    Optimize(mles, likelihoodFunction);

    results = estimators.ExtractMLEs( & likelihoodFunction, model_id_to_object);

	// Ensure frequencies were optimizes 
    assert(Abs(results["global"]) == 3, "frequencies not optimized");
	
}

test_binaryML()							 
							 
						 
												  