# Reinforcement-Learning

<h1>Optimizing collective fieldtaxis of swarming agents through reinforcement learnng</h1>
<h3>Introduction</h3><br>
The swarming of animal groups has fascinated scientists in fields ranging from
biology to physics to engineering. The complex swarming patterns often arise from
simple interactions between individuals to the benefit of the collective whole. Here we
show that a machine-learning technique can be employed to tune these underlying
parameters and optimize the resulting performance.
In this project we try to simulate a group of fish which try to find preferred regions
in a sea or ocean where there might be less light so that the predator doesn’t find them.
They do it collectively by interacting with each other on an individual level that benefits
the whole swarm.
<br>
<h3>Algorithm and Simulations</h3>
Part - I
<br><br>
Firstly I casted a static light field using the equation<br><br>
<h4>
  𝐹(𝑋) =𝑘∑[𝐴<sub>𝑘</sub>𝑐𝑜𝑠(𝑘. 𝑋) + 𝐵𝑘𝑠𝑖𝑛(𝑘. 𝑋)]
</h4>
<br>
We do the sum over wave vectors, k =(𝑥,𝑘𝑦) 𝑘
𝑘 which runs through
𝑥
,𝑘
𝑦
2π
𝐿
We get 𝐴 and using gaussian normal distribution with mean of 0.5 and standard 𝑘
𝐵𝑘
deviation 0.05
Then I initialized positions of the agents uniformly within a circle of radius R = ,
𝑁
π
where N = number of agents , here N = 16
In fishes we can see a coherent flocking happens by interacting through repulsion at
short range 𝑟 , orientation at intermediate range , and attraction at long range .
𝑟
𝑟
𝑜
𝑟
𝑎
With each time step ∆𝑡, we update the positions and velocities as
𝑥 =
𝑖
(𝑡 + ∆𝑡) 𝑥
𝑖
(𝑡) +
1
2
[𝑣
𝑖
(𝑡) + 𝑣
𝑖
(𝑡 + ∆𝑡)]
𝑣 = ]
𝑖
(𝑡 + ∆𝑡) [𝑣𝑚𝑎𝑥
𝐹(𝑥
𝑖
(𝑡))] 𝑅
𝑛𝑜𝑖𝑠𝑒
𝑑
𝑖
(𝑡)
The position of particles changes accordingly with the velocity vector
And the magnitude of velocity of particles changes with the amount light received
And the direction of the particles changes accordingly with the following conditions
i) If there are any neighbors within the zone of repulsion 𝑟 , they need to repel so
𝑟
ii) If there neighbors in the zone of repulsion but within zone of orientation 𝑟 ,then they
𝑜
need orient in the same direction. If within zone of attraction 𝑟 , then they need to
𝑎
attract so
iii) If there are no neighbors within zone of attraction 𝑟 , then there will be no change in
𝑎
the direction so
In training we allow the agents to move up to time t = 𝑡 where = 100 in our case *
𝑡
*
We can see the following path by the agents
Here we chose 𝑟 = 1.95 , = 2.05 & = 1.00
𝑜
𝑟
𝑎
𝑟
𝑟
We need to optimize 𝑟 using reinforcement learning so that the agents can find the
𝑜
, 𝑟
𝑎
darkest spot easily.
For this we run the training algorithm for 𝑁 times indexed by .
𝑡𝑟𝑎𝑖𝑛
α = 1, 2,.... 𝑁
𝑡𝑟𝑎𝑖𝑛
We take the average field intensity perceived by the agents so that we can use this as
reward
Then we define the reward as
Q = max{ 0 ,𝑓(0) - 𝑓(𝑡 }
*
)
So, if f(t) decreases over time 𝑡 then the reward is going to increase. *
At α-th training session with initial positions and velocities of agents, we evaluate the
rewards at three nearby points: 𝑄 at ( ), at and at
0
𝑟
𝑜
(α)
, 𝑟
𝑎
(α) 𝑄1
(𝑟
𝑜
(α) − δ, 𝑟
𝑎
(α)
) 𝑄2
(𝑟 with deviation . And we update the parameters as
𝑜
(α)
, 𝑟
𝑎
(α) + δ) δ = (α + 1)
−1/4
After many iterations we can see how 𝑟 are changing and get the optimised .
𝑜
, 𝑟
𝑎
𝑟
𝑜
, 𝑟
𝑎
We can see that 𝑟 are approaching 1.75, 1.68.
𝑜
, 𝑟
𝑎
Part - II
We need to train our agents with allowed continuum values of ro and ra i.e., [1,3] .
We need to train on each and every pair of ro and ra values and take the average
reward over different static fields.
After training over 1600 fields we got
In this we can see that ra greater than ro will give us more reward.
Conclusion
In this project I have simulated a static light field and introduced agents to swarm
around and find the darkest place. We can see how the simple interactions between
agents benefits the whole swarm. I have implemented an algorithm to optimize 𝑟
𝑜
, 𝑟
𝑎
using reinforcement learning which was given in the paper. We can see how optimizing
𝑟 can help in effectively sensing the gradient.
