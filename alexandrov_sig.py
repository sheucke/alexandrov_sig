from SigProfilerExtractor import sigpro as sig


def main_function():
    # to get input from vcf files
    #path_to_example_folder_containing_vcf_files = sig.importdata("vcf")
    # you can put the path to your folder containing the vcf     samples
    data = "./vcf_files"
    sig.sigProfilerExtractor("vcf", "alexandrov_results",
                             data, minimum_signatures=1, maximum_signatures=3)


if __name__ == "__main__":
    main_function()
