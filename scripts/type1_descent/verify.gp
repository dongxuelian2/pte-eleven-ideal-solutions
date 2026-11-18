default(parisizemax,1000000000);
read("scripts/type1_descent/selmer.gp");
read("scripts/type1_descent/tangents.gp");
read("scripts/type1_descent/pairing.gp");
verify(b,reps,points,expected)={
  my(es=rootsE(b),dim=selmerdim(es),M=pairmatrix(es,reps,points));
  if(dim!=expected || matrank(M)!=4,error("descent mismatch"));
  print("b0=",b,", Selmer dimension ",dim,", pairing minor ",lift(M),", rank <= ",dim-6);
};
verify(0,reps0,points0,8);
verify(1,reps1,points1,7);
