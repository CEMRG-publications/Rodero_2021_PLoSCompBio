1)
in "data" folder, put your folder with given foldername, e.g. "tag"

ls data/tag/

X.txt
Y.txt
xlabels.txt
ylabels.txt

2)
in "tmp" folder, mkdir a folder with same foldername as the one containing your data,
and inside this folder mkdir as many empty folders as the number of features to be emulated (corresponding to len(ylabels), counting from 0):

mkdir tmp/tag
cd tmp/tag
for i in {0..7}; do mkdir ${i}; done

(for 8 output features to be emulated in this case)

3)

./run_all_gpe.sh

NOTE: before running this, if you are not planning to use GPU, comment the line indicated in the file "utils/concurrent.py"

This script automatically runs all the GPEs training in parallel in a 5-fold cross-validation, your PC should manage to hold the computational stress. If it doesn't, you may have to reduce the number of cross-validation folds that start simultaneously (email me).

4)

When (3) is completed you can run GSA, this will use the best scoring GPEs from the cross-validation.

./run_all_gsa.sh

