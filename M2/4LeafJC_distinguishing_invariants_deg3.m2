-- 4-leaf JC 4-cycle networks
-- The ideal is computed for one labeling of the topology, with leaf 1 at the reticulation vertex

restart;

A = (1,1);C = (1,-1);G = (-1,1); T = (-1,-1);

n = hashTable{0 => A, 1 => C, 2 => G, 3 => T}
x = hashTable{A => 0, C => 1, G => 1, T => 1}--JC
y = hashTable{A => 0, C => 1, G => 2, T => 3}

L = {};
for i from 0 to 3 do (
      for j from 0 to 3 do (
            for k from 0 to 3 do (
                  for l from 0 to 3 do (
                        -- only consistent leaf labellings!
                        if n#i#0 * n#j#0 * n#k#0 * n#l#0 == 1 and n#i#1 * n#j#1 * n#k#1 * n#l#1 == 1 then L = append(L, (n#i, n#j, n#k, n#l));
                  )
            )
      )
)

L1 = toList apply(L,i->q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)));
L2 = {b_0,b_1,c_0,c_1,d_0,d_1,e_0,e_1,
      f_0,f_1,g_0,g_1,h_0,h_1};

R = QQ[(L2|L1), Degrees => toList(#L2 : 1) | toList(#L1 : 6), MonomialOrder => {#L2, #L1}];

--4-cycle Network Topology
Eqn0123 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#1))*c_(x#(i#2))*d_(x#(i#3))*f_(x#(i#0))*g_(x#(i#0#0*i#1#0,i#0#1*i#1#1))*h_(x#(i#0#0*i#1#0*i#2#0,i#0#1*i#1#1*i#2#1)) - 
b_(x#(i#1))*c_(x#(i#2))*d_(x#(i#3))*e_(x#(i#0))*g_(x#(i#1))*h_(x#(i#1#0*i#2#0,i#1#1*i#2#1)));
      
I0123 = ideal(Eqn0123);
G = selectInSubring(1,gens gb(I0123, DegreeLimit => 18));
Gl = flatten entries G;
Gm = new MutableList from Gl;
counter = new MutableList from apply(Gl, f -> 1);

Eqn0132 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#1))*c_(x#(i#3))*d_(x#(i#2))*f_(x#(i#0))*g_(x#(i#0#0*i#1#0,i#0#1*i#1#1))*h_(x#(i#0#0*i#1#0*i#3#0,i#0#1*i#1#1*i#3#1)) - 
b_(x#(i#1))*c_(x#(i#3))*d_(x#(i#2))*e_(x#(i#0))*g_(x#(i#1))*h_(x#(i#3#0*i#1#0,i#3#1*i#1#1)));

I0132 = ideal(Eqn0132);
G1 = apply(Gl, f -> f % I0132);
for i from 0 to #G1-1 do if G1_i == 0 then counter#i = counter#i + 1;

Eqn0213 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#2))*c_(x#(i#1))*d_(x#(i#3))*f_(x#(i#0))*g_(x#(i#0#0*i#2#0,i#0#1*i#2#1))*h_(x#(i#0#0*i#2#0*i#1#0,i#0#1*i#2#1*i#1#1)) - 
b_(x#(i#2))*c_(x#(i#1))*d_(x#(i#3))*e_(x#(i#0))*g_(x#(i#2))*h_(x#(i#2#0*i#1#0,i#2#1*i#1#1)));

I0213 = ideal(Eqn0213)
G2 = apply(Gl, f -> f % I0213)
for i from 0 to #G2-1 do if G2_i == 0 then counter#i = counter#i + 1;

Eqn1023 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#0))*c_(x#(i#2))*d_(x#(i#3))*f_(x#(i#1))*g_(x#(i#0#0*i#1#0,i#0#1*i#1#1))*h_(x#(i#0#0*i#1#0*i#2#0,i#0#1*i#1#1*i#2#1)) - 
b_(x#(i#0))*c_(x#(i#2))*d_(x#(i#3))*e_(x#(i#1))*g_(x#(i#0))*h_(x#(i#0#0*i#2#0,i#0#1*i#2#1)));

I1023 = ideal(Eqn1023)
G3 = apply(Gl, f -> f % I1023)
for i from 0 to #G3-1 do if G3_i == 0 then counter#i = counter#i + 1;

Eqn1032 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#0))*c_(x#(i#3))*d_(x#(i#2))*f_(x#(i#1))*g_(x#(i#0#0*i#1#0,i#0#1*i#1#1))*h_(x#(i#0#0*i#1#0*i#3#0,i#0#1*i#1#1*i#3#1)) - 
b_(x#(i#0))*c_(x#(i#3))*d_(x#(i#2))*e_(x#(i#1))*g_(x#(i#0))*h_(x#(i#0#0*i#3#0,i#0#1*i#3#1)));

I1032 = ideal(Eqn1032)
G4 = apply(Gl, f -> f % I1032)
for i from 0 to #G4-1 do if G4_i == 0 then counter#i = counter#i + 1;

Eqn1203 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#2))*c_(x#(i#0))*d_(x#(i#3))*f_(x#(i#1))*g_(x#(i#1#0*i#2#0,i#1#1*i#2#1))*h_(x#(i#1#0*i#2#0*i#0#0,i#1#1*i#2#1*i#0#1)) - 
b_(x#(i#2))*c_(x#(i#0))*d_(x#(i#3))*e_(x#(i#1))*g_(x#(i#2))*h_(x#(i#2#0*i#0#0,i#2#1*i#0#1)));

I1203 = ideal(Eqn1203)
G5 = apply(Gl, f -> f % I1203)
for i from 0 to #G5-1 do if G5_i == 0 then counter#i = counter#i + 1;

Eqn2013 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#0))*c_(x#(i#1))*d_(x#(i#3))*f_(x#(i#2))*g_(x#(i#2#0*i#0#0,i#2#1*i#0#1))*h_(x#(i#2#0*i#0#0*i#1#0,i#2#1*i#0#1*i#1#1)) - 
b_(x#(i#0))*c_(x#(i#1))*d_(x#(i#3))*e_(x#(i#2))*g_(x#(i#0))*h_(x#(i#0#0*i#1#0,i#0#1*i#1#1)));

I2013 = ideal(Eqn2013)
G6 = apply(Gl, f -> f % I2013)
for i from 0 to #G6-1 do if G6_i == 0 then counter#i = counter#i + 1;

Eqn2031 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#0))*c_(x#(i#3))*d_(x#(i#1))*f_(x#(i#2))*g_(x#(i#2#0*i#0#0,i#2#1*i#0#1))*h_(x#(i#2#0*i#0#0*i#3#0,i#2#1*i#0#1*i#3#1)) - 
b_(x#(i#0))*c_(x#(i#3))*d_(x#(i#1))*e_(x#(i#2))*g_(x#(i#0))*h_(x#(i#0#0*i#3#0,i#0#1*i#3#1)));

I2031 = ideal(Eqn2031)
G7 = apply(Gl, f -> f % I2031)
for i from 0 to #G7-1 do if G7_i == 0 then counter#i = counter#i + 1;

Eqn2103 = toList apply(L, i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#1))*c_(x#(i#0))*d_(x#(i#3))*f_(x#(i#2))*g_(x#(i#2#0*i#1#0,i#2#1*i#1#1))*h_(x#(i#2#0*i#1#0*i#0#0,i#2#1*i#1#1*i#0#1)) - 
b_(x#(i#1))*c_(x#(i#0))*d_(x#(i#3))*e_(x#(i#2))*g_(x#(i#1))*h_(x#(i#1#0*i#0#0,i#1#1*i#0#1)));

I2103 = ideal(Eqn2103)
G8 = apply(Gl, f -> f % I2103)
for i from 0 to #G8-1 do if G8_i == 0 then counter#i = counter#i + 1;

Eqn3012 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#0))*c_(x#(i#1))*d_(x#(i#2))*f_(x#(i#3))*g_(x#(i#3#0*i#0#0,i#3#1*i#0#1))*h_(x#(i#3#0*i#0#0*i#1#0,i#3#1*i#0#1*i#1#1)) - 
b_(x#(i#0))*c_(x#(i#1))*d_(x#(i#2))*e_(x#(i#3))*g_(x#(i#0))*h_(x#(i#0#0*i#1#0,i#0#1*i#1#1)));

I3012 = ideal(Eqn3012)
G9 = apply(Gl, f -> f % I3012)
for i from 0 to #G9-1 do if G9_i == 0 then counter#i = counter#i + 1;

Eqn3021 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#0))*c_(x#(i#2))*d_(x#(i#1))*f_(x#(i#3))*g_(x#(i#3#0*i#0#0,i#3#1*i#0#1))*h_(x#(i#3#0*i#0#0*i#2#0,i#3#1*i#0#1*i#2#1)) - 
b_(x#(i#0))*c_(x#(i#2))*d_(x#(i#1))*e_(x#(i#3))*g_(x#(i#0))*h_(x#(i#0#0*i#2#0,i#0#1*i#2#1)));

I3021 = ideal(Eqn3021)
G10 = apply(Gl, f -> f % I3021)
for i from 0 to #G10-1 do if G10_i == 0 then counter#i = counter#i + 1;

Eqn3102 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
b_(x#(i#1))*c_(x#(i#0))*d_(x#(i#2))*f_(x#(i#3))*g_(x#(i#3#0*i#1#0,i#3#1*i#1#1))*h_(x#(i#3#0*i#1#0*i#0#0,i#3#1*i#1#1*i#0#1)) - 
b_(x#(i#1))*c_(x#(i#0))*d_(x#(i#2))*e_(x#(i#3))*g_(x#(i#1))*h_(x#(i#1#0*i#0#0,i#1#1*i#0#1)));

I3102 = ideal(Eqn3102)
G11 = apply(Gl, f -> f % I3102)
for i from 0 to #G11-1 do if G11_i == 0 then counter#i = counter#i + 1;

for i from 0 to #counter - 1 do if counter#i == 12 then Gm#i = 0; 

phyInv = toList(Gm)
"4LeafJC_GB_distinguishing.txt" << toString(phyInv) << close;
quit;

