awk 'BEGIN{
natoms=12
for(i=0;i<10;i++) {
print natoms
print "0 0 0"
for(j=0;j<natoms;j++) {
  print("X",i+j*j+1,0.0,0.0)
}
}
}'
