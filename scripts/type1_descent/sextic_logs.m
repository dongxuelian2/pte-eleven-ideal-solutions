// Read the same abscissas as the PARI four-cover verification.
SetColumns(0);
Q:=Rationals(); Z:=Integers(); A<t>:=PolynomialRing(Q);
K<w>:=NumberField(t^6-2*t^5-171*t^4+7227*t^3-60030*t^2+184672*t-195552);
r:=7925/28493184*w^5+13037/35616480*w^4-1332239/28493184*w^3
   +263051489/142465920*w^2-190118627/17808240*w+47832223/2968040;
h:=-800*t^6+2625*t^5-12050*t^4+529*t^3+54448*t^2+53744*t+14400;
assert Evaluate(h,r) eq 0;
data:=eval Read("scripts/type1_descent/sextic_local_points.txt");

function Coordinates(z,f)
    a:=Eltseq(z);
    if f eq 3 then a:=Eltseq(a[1]);end if;
    assert forall{t:t in a|AbsolutePrecision(t) ge 24};
    return [Z!v mod 2^24:v in a];
end function;

for item in Decomposition(MaximalOrder(K),2) do
    p:=item[1];e:=item[2];f:=InertiaDegree(p);
    if e*f eq 1 then continue;end if;
    i:=e eq 2 select 1 else 2;
    L,em:=Completion(K,p:Precision:=500);
    c:=2^(-Valuation(em(r)) div e);
    roots:=[K|0,-5*r^2*c^2,4*r*c^2];
    E:=EllipticCurve([L|0,-em(&+roots),0,em(roots[2]*roots[3]),0]);
    zs:=[];
    for row in data[i] do
        x:=K!row;
        if x in roots then
            Append(~zs,[0:k in [1..e*f]]);
        else
            u:=em(x);ok,v:=IsSquare(&*[u-em(a):a in roots]);assert ok;
            P:=E![u,v,1];W:=(105*2^7)*P;assert W ne E!0;
            z:=-W[1]/W[2];assert Valuation(z) ge 7*e;
            Append(~zs,Coordinates(z,f));
        end if;
    end for;
    bas:=i eq 1 select [2,3] else [3,4,5];
    B:=Matrix(Q,[zs[j]:j in bas]);Binv:=B^-1;
    assert forall{v:v in Eltseq(256*Binv)|Valuation(v,2) ge 0};
    D:=Matrix(Q,zs[1..#zs-1])*Binv;
    ref:=Vector(Q,zs[#zs])*Binv/2;
    ray:=256*Vector(Q,Coordinates(em(1/(c*r)),f))*Binv;
    assert forall{v:v in Eltseq(D) cat Eltseq(ref) cat Eltseq(ray)|Valuation(v,2) ge 0};
    F:=GF(2);
    expected:=i eq 1 select [[0,1,0,0],[0,0,1,0]]
                      else [[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]];
    assert ChangeRing(Transpose(D),F) eq Matrix(F,expected);
    assert Eltseq(ChangeRing(ref,F)) eq (i eq 1 select [1,0] else [1,1,1]);
    assert Eltseq(ChangeRing(ray,F)) eq [1:j in [1..e*f]];
    print "Verified dyadic logarithms for ramification degree",e,"and residue degree",f;
end for;
