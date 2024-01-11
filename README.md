# Stress-based-continuous-planar-path-planning-for-additive-manufacturing

**Paper:** 
[https://doi.org/10.1016/j.advengsoft.2023.103544]

**Abstract:**
This paper presents a novel continuous planar path planning algorithm based on principal stress orientation field to improve structural performance of the products fabricated with additive manufacturing. For a given 3D shape, it is first sliced into planar layers with specified layer thickness, which are discretized into grids for the optimization of tool-paths. Then, the principal stresses and the corresponding orientation field are calculated with the finite element analysis method, which are projected onto grids. Since the principal stress orientations of elements with principal stress values close to zero are arbitrary, a weighted filtering strategy is applied to generate a smooth principal stress orientation field. Finally, the tool-paths are generated with genetic algorithm through minimizing the idle travel distance of the nozzle. Furthermore, the tool-path width is adaptively adjusted according to the tool-path orientation to avoid the overfill issue caused by the discretization. Additionally, a method is introduced to achieve global continuous tool-paths, aiming at eliminating idle travel completely. The optimized tool-paths are simulated and compared with the state-of-the-art methods. Experimental results demonstrate that the objects fabricated with the proposed method have superior structural performance.

**Contact information:**
Hao Xu (haoxu5640@gmail.com)
