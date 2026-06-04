mkdir -p $PREFIX/bin

cp $SRC_DIR/gfa_stats.py $PREFIX/bin/
cp $SRC_DIR/gfa_to_fasta.py $PREFIX/bin/
cp $SRC_DIR/lrSAGA.py $PREFIX/bin/
cp $SRC_DIR/preprocess_mbg_gfa.py $PREFIX/bin/

chmod +x $PREFIX/bin/gfa_stats.py
chmod +x $PREFIX/bin/gfa_to_fasta.py
chmod +x $PREFIX/bin/lrSAGA.py
chmod +x $PREFIX/bin/preprocess_mbg_gfa.py

# lrSAGA utils
mkdir -p $SP_DIR/utils
cp $SRC_DIR/utils/__init__.py $SP_DIR/utils/
cp $SRC_DIR/utils/*.py $SP_DIR/utils/
