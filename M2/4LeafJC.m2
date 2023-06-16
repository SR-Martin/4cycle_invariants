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
L2 = {a_0,a_1,b_0,b_1,c_0,c_1,d_0,d_1,e_0,e_1,
      f_0,f_1,g_0,g_1,h_0,h_1};

R = QQ[(L2|L1), MonomialOrder=> Eliminate 16];

--4-cycle Network Topology
Eqn_4 = toList apply(L,i->
q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)) - 
a_(x#(i#0))*b_(x#(i#1))*c_(x#(i#2))*d_(x#(i#3))*f_(x#(i#0))*g_(x#(i#0#0*i#1#0,i#0#1*i#1#1))*h_(x#(i#3)) - 
a_(x#(i#0))*b_(x#(i#1))*c_(x#(i#2))*d_(x#(i#3))*e_(x#(i#0))*g_(x#(i#1))*h_(x#(i#0#0*i#3#0,i#0#1*i#3#1)));
      

K_4 = ideal(Eqn_4);

"4LeafJC_GB.txt" << selectInSubring(1,gens gb(K_4)) << close;
quit;
