# setting max db size to 20Gb. kraken tries to load the whole database into RAM and crashes with a memmory allocation error
# in this version. Rather than rolling back versions, the devs reccomend just setting a max db size limit. Kraken can then
# process in chunks that fit into memory.
kraken2-build --max-db-size 20000000000 --standard --db kraken_db
