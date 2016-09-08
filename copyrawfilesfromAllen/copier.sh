#!/usr/bin/bash
listoflibs="1146 1147 1148 1149 1131 1132 1163 
	1164 1178 1182 1183 1184 1186 1188 1190 1191 
	1192 1150 1162 1165 1166 1167 1168 1169"
path="/var/www/subdirectories_for_interface/temp_output/library_"

#copy
for number in $listoflibs
do 
	newpath=$path$number
        cp -rf $newpath TRANSFER/
done

#zipped
tar czf transfer.tgz TRANSFER

#transfer/copied
scp -P 1657 transfer.tgz modupeore17@geco.iplantc.org:~/modupeore17/
