sed "s/EXPNUMBER/$1/" input_norestart_template >/tmp/startzoupa
./zoupa </tmp/startzoupa >$2
