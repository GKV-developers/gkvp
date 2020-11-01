#!/bin/sh

old=f0.51
new=f0.52

rm -f rename.log
touch rename.log

echo "### Replace *${old}* to *${new}* ###" >> rename.log
grep -r ${old} ./extra_tools/* >> rename.log
grep -r ${old} ./lib/*         >> rename.log
grep -r ${old} ./run/*         >> rename.log
grep -r ${old} ./src/*         >> rename.log

  grep -rl ${old} ./extra_tools/* | xargs sed -i "s/${old}/${new}/g"
  grep -rl ${old} ./lib/*         | xargs sed -i "s/${old}/${new}/g"
  grep -rl ${old} ./run/*         | xargs sed -i "s/${old}/${new}/g"
  grep -rl ${old} ./src/*         | xargs sed -i "s/${old}/${new}/g"


echo "### Rename files from *${old}* to *${new}* ###" >> rename.log
find ./extra_tools/* -name *${old}* >> rename.log
find ./lib/*         -name *${old}* >> rename.log
find ./run/*         -name *${old}* >> rename.log
find ./src/*         -name *${old}* >> rename.log
echo >> rename.log

  find ./extra_tools/* -name *${old}* | xargs rename ${old} ${new}
  find ./lib/*         -name *${old}* | xargs rename ${old} ${new}
  find ./run/*         -name *${old}* | xargs rename ${old} ${new}
  find ./src/*         -name *${old}* | xargs rename ${old} ${new}


