import os
import fnmatch
import shutil

##########
# you can run this script after you downloaded new data:
# 
# the folder directory structure should  be as follows:
# 20240101_DARLIN_Mon_Gr_XXX
#     |
#     |___ Data
#         |___Gr_A1
#             |___Gr_A1_R1.fq.gz
#             |___Gr_A1_R2.fq.gz
# or
# 20240101_DARLIN_Mon_Gr_XXX
#     |
#     |___ raw_fastq
#         |___Gr_A1
#             |___Gr_A1_R1.fq.gz
#             |___Gr_A1_R2.fq.gz

# run rename script in the folder 20240101_DARLIN_Mon_Gr_XXX:
# cd 20240101_DARLIN_Mon_Gr_XXX then run this rename script
##########


# get files
def find_files(directory, pattern):
    matched_files = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, pattern):
            matched_files.append(os.path.join(root, filename))
    matched_files = [item for item in matched_files if 'checkpoint' not in item]
    return matched_files

def process_format(pattern, remove_str, replace_str, process_folders):
    files = find_files(cwd, pattern)

    if len(files):
        for item in files:
            file_dir = os.path.dirname(item)
            file_basename = os.path.basename(item)
            filter_name = file_basename.replace(remove_str, replace_str)
            shutil.move(item, os.path.join(file_dir,filter_name))

        if process_folders:
            folders = set(os.path.dirname(item) for item in files)
            for item in folders:
                filter_name = item.replace('-', '_')
                shutil.move(item, filter_name)

    # rename folder: Data -> raw_astq
    if os.path.exists(os.path.join(cwd, 'Data')):
        shutil.move(os.path.join(cwd, 'Data'), os.path.join(cwd, 'raw_fastq'))

def main(cwd):
    patterns = ['*.fq.gz', '*.fastq.gz', '*.fastq.gz', '*.fastq.gz']
    remove_strs = ['fq.gz', '-', '_L001', '_001']
    replace_strs = ['fastq.gz', '_', '', '']


    for remove_str in remove_strs:
        process_folders = False
        index = remove_strs.index(remove_str)
        pattern = patterns[index]
        replace_str = replace_strs[index]
        if index == len(patterns)-1:
            process_folders = True
        process_format(pattern, remove_str, replace_str, process_folders)

if __name__ == '__main__':
    cwd = os.getcwd()
    main(cwd)
