
echo "-------Basic relative--------------------------"
../src/fastdate --method_relative --tree_file ../example/small.tre  --out_file rel_dates.tre --bd_mu 1 --bd_lambda 4  --bd_rho 0.01 --show_tree --grid 100 --rate_mean 5 --rate_variance 1 --bd_psi 0.1

echo "-------Basic node prior------------------------"
../src/fastdate --method_nodeprior --tree_file ../example/small.tre --prior_file ../example/node_prior.txt --out_file node_prior.tre --bd_mu 1 --bd_lambda 4 --max_age 100 --bd_rho 1 --show_tree --grid 100 --rate_mean 5 --rate_variance 1

echo "-------Basic tip dates-------------------------"
../src/fastdate --method_tipdates --tree_file ../example/small.tre --prior_file ../example/tip_prior.txt --out_file tip_dates.tre --bd_mu 1 --bd_lambda 4 --max_age 100 --bd_rho 0.01 --show_tree --grid 100 --rate_mean 5 --rate_variance 1 --bd_psi 0.1

echo "-------Sample relative--------------------------"
../src/fastdate --method_relative --tree_file ../example/small.tre  --out_file rel_dates.tre --bd_mu 1 --bd_lambda 4  --bd_rho 0.01 --show_tree --grid 100 --rate_mean 5 --rate_variance 1 --bd_psi 0.1 --sample 1
if [ $(diff rel_dates.tre rel_dates.tre.sampled | wc -w) -eq 0 ];
then 
  echo "Sample doesn't differ from best tree"
fi

echo "-------Sample node prior-------------------------"
../src/fastdate --method_nodeprior --tree_file ../example/small.tre --prior_file ../example/node_prior.txt --out_file node_prior.tre --bd_mu 1 --bd_lambda 4 --max_age 100 --bd_rho 1 --show_tree --grid 100 --rate_mean 5 --rate_variance 1 --sample 1

if [ $(diff node_prior.tre node_prior.tre.sampled | wc -w) -eq 0 ];
then 
  echo "Sample doesn't differ from best tree"
fi


echo "--------Sample tip dates-------------------------"
../src/fastdate --method_tipdates --tree_file ../example/small.tre --prior_file ../example/tip_prior.txt --out_file tip_dates.tre --bd_mu 1 --bd_lambda 4 --max_age 100 --bd_rho 0.01 --show_tree --grid 100 --rate_mean 5 --rate_variance 1 --bd_psi 0.1 --sample 1

if [ $(diff tip_dates.tre tip_dates.tre.sampled | wc -w) -eq 0 ];
then 
  echo "Sample doesn't differ from best tree"
fi


echo "-------Interval relative--------------------------"
../src/fastdate --method_relative --tree_file ../example/small.tre  --out_file rel_dates.tre --bd_mu 1 --bd_lambda 4  --bd_rho 0.01 --show_tree --grid 100 --rate_mean 5 --rate_variance 1 --bd_psi 0.1 --cred_interval 0.95

echo "-------Interval node prior------------------------"
../src/fastdate --method_nodeprior --tree_file ../example/small.tre --prior_file ../example/node_prior.txt --out_file node_prior.tre --bd_mu 1 --bd_lambda 4 --max_age 100 --bd_rho 1 --show_tree --grid 100 --rate_mean 5 --rate_variance 1 --cred_interval 0.95

echo "-------INterval tip dates-------------------------"
../src/fastdate --method_tipdates --tree_file ../example/small.tre --prior_file ../example/tip_prior.txt --out_file tip_dates.tre --bd_mu 1 --bd_lambda 4 --max_age 100 --bd_rho 0.01 --show_tree --grid 100 --rate_mean 5 --rate_variance 1 --bd_psi 0.1 --cred_interval 0.95

