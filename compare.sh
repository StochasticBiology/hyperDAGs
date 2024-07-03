g++ -o3 hyperhmm-friends.cpp -o hyperhmm-friends.ce -g
./hyperhmm-friends.ce --obs test.txt --longitudinal --nboot 0
./hyperhmm-friends.ce --obs test.txt --longitudinal --nboot 0 --edgelist test-table.txt

head -n100 tb_drug.txt > tb-test.txt
./hyperhmm-friends.ce --obs tb-test.txt --longitudinal --eps 1e-3 --nboot 0 --label tb-normal
./hyperhmm-friends.ce --obs tb-test.txt --longitudinal --eps 1e-3 --nboot 0 --edgelist TB-table-ss.txt --label tb-constrained
./hyperhmm-friends.ce --obs tb-test.txt --longitudinal --eps 1e3 --nboot 0 --label tb-normal-one
./hyperhmm-friends.ce --obs tb-test.txt --longitudinal --eps 1e3 --nboot 0 --edgelist TB-table-ss.txt --label tb-constrained-one

./hyperhmm-friends.ce --obs tb_drug.txt --longitudinal --eps 1e-3 --nboot 0 --label tb-full-normal
./hyperhmm-friends.ce --obs tb_drug.txt --longitudinal --eps 1e3 --nboot 0 --edgelist TB-table-ss.txt --label tb-full-constrained-one

% TB-table.txt comes from the hyperDAGs R code

