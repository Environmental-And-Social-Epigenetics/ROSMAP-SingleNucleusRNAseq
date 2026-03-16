source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh

cd /om/scratch/Mon/mabdel03/SocialIsolation/

#python ver = 3.9

#python3 -m pip install git+https://github.com/yoseflab/Compass.git --upgrade

conda activate /om2/user/mabdel03/conda_envs/cplex_env2

export CPLEX_STUDIO_DIR=/om/scratch/Mon/mabdel03/SocialIsolation/opt/ibm/ILOG/CPLEX_Studio2211
export PATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PATH
export PYTHONPATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PYTHONPATH

echo "y" | compass --data dejpseudobulk_male_Mic.tsv --num-processes 10 --species homo_sapiens --output-dir CompassPMMic

echo "y" | compass --data dejpseudobulk_female_Mic.tsv --num-processes 10 --species homo_sapiens --output-dir CompassPFMic