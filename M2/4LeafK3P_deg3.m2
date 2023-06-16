-- 4-leaf K3P 4-cycle networks
-- The ideal is computed for one labeling of the topology, with leaf 0 at the reticulation vertex

restart;

A = (1,1);C = (1,-1);G = (-1,1); T = (-1,-1);

n = hashTable{0 => A, 1 => C, 2 => G, 3 => T}
x = hashTable{A => 0, C => 1, G => 2, T => 3}--K3P
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
L2 = {b_1,b_2,b_3,c_1,c_2,c_3,d_1,d_2,d_3,e_1,e_2,e_3,f_1,f_2,f_3,g_1,g_2,g_3,h_1,h_2,h_3};

R = QQ[({t} | L2 | L1), Degrees => toList(#L2 + 1 : 1) | toList(#L1 : 6), MonomialOrder => {#L2 + 1, #L1}];

--4-cycle Network Topology
Eqn_4 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - (if x#(i#1) == 0 then t else b_(x#(i#1))) * (if x#(i#2) == 0 then t else c_(x#(i#2))) * (if x#(i#3) == 0 then t else d_(x#(i#3))) * ( (if x#(i#0) == 0 then t else f_(x#(i#0))) * (if x#(i#0#0*i#1#0,i#0#1*i#1#1) == 0 then t else g_(x#(i#0#0*i#1#0,i#0#1*i#1#1))) * (if x#(i#3) == 0 then t else h_(x#(i#3))) + (if x#(i#0) == 0 then t else e_(x#(i#0))) * (if x#(i#1) == 0 then t else g_(x#(i#1))) * (if x#(i#0#0*i#3#0,i#0#1*i#3#1) == 0 then t else h_(x#(i#0#0*i#3#0,i#0#1*i#3#1)))));

K_4 = ideal(Eqn_4);
G = selectInSubring(1,gens gb(K_4, DegreeLimit => 18));
"4LeafK3P_GB_deg3.txt" << G << close;

quit;
