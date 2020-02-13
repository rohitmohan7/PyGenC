# PyGenC

![istockphoto-935317640-612x612](https://user-images.githubusercontent.com/40918560/74388306-16880680-4dc9-11ea-96d2-244801f38682.jpg)

Are you utilizing the full functionality of the C/C++ engine efficiently to code? Are you following true test driven Development (TDD) when you are writing both the unit tests and the code?........Introducing PygenC! The world’s first repo that generates Code automatically from Unit test definitions using Darwin’s mechanism of evolution:  

The basic idea of PygenC is to generate a base repo from a test definition file with many unknowns such as function execution, inputs, return type etc. PygenC parses the entire gcc/g++ library to see the available C functions or ‘genetic mutations’ that can be introduced to the ‘genes’ or repo without causing a compile error. An initial population is generated with random seed of these functions up to the specified carrying capacity of the population and natural/artificial ‘selection’ is performed. The mutants in the population are given a fitness score based on the Unit tests that pass and every generation half of the mutants with the lowest fitness score perishes. Of the mutants that survive, based on the type of reproduction selected, binary fission (simple, fast) or sex (complex), mutants either split into two mutants/repos with mutation, or produce gametes with exactly half of the genetic information/code with mutation and crossover happening between these gametes. The gametes of the best fit mutants are selected for mating to produce the next fit mutant until the best fit mutant is found. PygenC can be used to evolve existing repos.

NOTE: Pygen is still under development

"Code should not be written, it should be evolved!"
