

stataMP -e do 1-format-neear-epi.do
stataMP -e do 2-format-adm-epi.do
stataMP -e do 3-format-mb-epi.do
stataMP -e do 4-append-epi-data.do
stataMP -e do 5-format-neear-wq.do
stataMP -e do 6-format-adm-wq.do
stataMP -e do 7-format-mb-wq.do
stataMP -e do 8-append-wq-data.do
stataMP -e do 9-avg-wq-data.do
stataMP -e do 10-make-analysis-dataset.do

