# toy

This is the AGB project by [**Álvaro Abella** :person_with_blond_hair:](https://github.com/alvaroabascar) and [**Clàudia Fontserè** :no_good:](https://github.com/claudefa)!

- [Organization of the code](#codeorganization)
    * [HMM.pm](#hmmpm)
    * [MODELS.pm](#modelspm)
    * [AUX.pm](#auxpm)
    * [countfreq.pl](#countfreqpl)
    * [test.pl](#testpl)
- [Objectives and solutions](#objectives-and-solutions)
    1. [Implementing a program to sample sequences.](#objective1)
    2. [Implementing Viterbi and testing with a toy model.](#objective2)
    3. [Consider the binding of U1 snRNP](#objective3)
    4. [Incorporate a TIA-1 binding site](#objective4)
    5. [Assessment of the performance](#objective5)

##Organization of the code
This project is split in different modules and scripts. Here is a little reference of what each file contains. An index
of functions from each module (devised for personal use, but maybe useful) can be found in [functions.pm](https://github.com/alvaroabascar/toy/blob/master/functions.pm).

#####[HMM.pm](https://github.com/alvaroabascar/toy/blob/master/HMM.pm)
Core functions needed to run a Hidden Markov Model, including the Viterbi algorithm, a function to run tests and some
abstractions to simplify the rest of the code.

#####[MODELS.pm](https://github.com/alvaroabascar/toy/blob/master/MODELS.pm)
The definition of each HMM model used along this project. Each model is composed by a set of states,
symbols, transition and emission probabilities.
- **toy**: toy model (with only one donor state)
- **U1**: adding most informative positions around the donor site
- **U1_all**: adding all positions around the donor site)
- **TIA1**: taking U1_all positions but adding one splice regulator
- **TIA2**: taking U1_all positions but adding two splice regulators
- **TIA3**: taking U1_all positions but adding one splice regulator with a different approach

You can find the model structure for them in the directory ./graphs

#####[AUX.pm](https://github.com/alvaroabascar/toy/blob/master/AUX.pm)
Contains some helper functions, not related to Markov models.

#####[countfreq.pl](https://github.com/alvaroabascar/toy/blob/master/countfreq.pl)
A simple script that we have used to calculate:

1. The frequencies of nucleotides at each position of the real ([real_human_donors](https://github.com/alvaroabascar/toy/blob/master/real_human_donors)) and false ([false_human_donors](https://github.com/alvaroabascar/toy/blob/master/false_human_donors)) donors.
2. From the frequencies in each position of both sets of sequences, we have calculated the Kullback-Leibler divergence among the distribution of nucleotides in the case of real and false donors, for each position.
3. From the Kullback-Leibler divergence we have calculated the *Jensen-Shannon divergence*, whose graph is shown [here](#objective3)

#####[sampling.pl](https://github.com/alvaroabascar/toy/blob/master/sampling.pl)
The functions needed to perform a sampling of *n* nucleotides from a given HMM.

#####[test.pl](https://github.com/alvaroabascar/toy/blob/master/sampling.pl)
A simple script to run a test, calling the corresponding function in HMM.pm


###Objectives and solutions
<a name="objective1"></a>
1. **Implement a program to sample sequences from an HMM**

We provide a simple script [sampling.pl](#samplingpl) which can be executed to produce a sequence of the desired length. Eg. to produce a sequence of 100 nucleotides:
```bash
./sampling.pl 100
```
**Note 1:** *the nucleotides are emitted according to the transition and emission probabilities of the HMM toy model. However, the sampling process will ignore the transition from intron to end in order to provide a sufficiently large sequence. If the requested length is very long, probably most of it will be intron.*

**Note 2:** *the default model to sample is the toy model. However, this can be altered by modifying the end of the script, specifically the line `HMM::select_model("toy")`. Eg. to sample from the HMM including the TIA-1 binding site:*
```perl
HMM::select_model("TIA");
print sampling(shift @ARGV), "\n";
```

<a name="objective2"></a>
2. **Implement the Viterbi algorithm and test it with the toy donor splice site model.**

The Viterbi algorithm is implemented as part of the module [HMM.pm](#hmmpm). The function is `HMM::viterbi`, which takes a sequence, builds a matrix of probabilities and pointers and returns the most probable path of states (for this task, it relies on the function `HMM::getPath`).

The testing is provided by the function `HMM::test`, which can be called from the script the script [test.pl](#testpl). This script can be easily modified to select a given model and a file of test sequences. The model is selected as shown in the [previous section](#objective1). The test is run calling `HMM::test('filename', n)`, where filename is a file containing a set of sequences, one per line, and n is the position (base 1) of this sequences where the donor site is found.

<a name="objective3"></a>
3. **Change the model into a donor site (5' splice site) model that considers the binding of the U1 snRNP, by extending the number of states that describe the exon-intron boundary. How many positions should you use for your model? Provide an argument for your answer.**

In order to find which positions we should use we calculated the *Jensen-Shannon* divergence among the sequences in [false_human_donors](https://github.com/alvaroabascar/toy/master/false_human_donors) and [real_human_donors](https://github.com/alvaroabascar/toy/master/real_human_donors), for each of the positions.

![Jensen-Shannon divergence](https://cloud.githubusercontent.com/assets/7307772/6461654/35be02d8-c1a2-11e4-820c-ebf03cc8441e.png)

From the previous graph it's clear that positions 3, 6 and 9 exhibit a large difference in the nucleotide frequencies. The positions 1, 2, 7 and 9 show a much smaller divergence, and positions 4 and 5 show no divergence at all. This last positions correspond to the dinucleotide **GT**, which is always seen in the donor site. The lack of divergence in this case is due to the fact that the [false_human_donors](https://github.com/alvaroabascar/toy/master/false_human_donors)) have been selected to have this dinucleotide in position 4-5. However, it is obvious that this two positions must be included in the model, as they are the most relevant.

Given this data we decided to take positions 3, 4, 5, 6, 7 and 8 (the 6th in order to be able to take the 8th). This model is included in [MODELS.pm](#modelspm) with the name "U1".

![U1 with the most relevant positions](https://cloud.githubusercontent.com/assets/7307772/6461646/1bd89c02-c1a2-11e4-9192-b7808d54e092.png)

We also decided to create another model taking into account all of the positions, which is called "U1_all".

![U1 with the all the positions](https://cloud.githubusercontent.com/assets/7307772/6461648/1eb2bed0-c1a2-11e4-880d-eff7ed907f6d.png)

<a name=runtestexample></a>
To test, for example, the first model, edit [test.pl](#testpl) as follows and run it:
```perl
HMM::select_model("U1");
HMM::test("testset_full.txt", 51);
```
This two models have provided, respectively, an accuracy of 86.58% and 94.00%.

<a name="objective4"></a>
4. **Incorporate into the previous model a state describing the presence of a TIA-1 binding site (a Uridine-rich sequence) immediately downstream of the donor site.**

We have incorporated this model under [MODELS.pm](#modelspm) with the name "TIA". After several tests we have seen that a good approach is to allow a small transition probability from donor to TIA, and a bigger probability of transition from donor to intron. Intron allows the transition both to itself and back to TIA, and TIA can lead to itself or to intron. This has provided a sensibility of around 91%, worse than the sensibility provided by "U1_all".

![TIA1](https://cloud.githubusercontent.com/assets/7307772/6461636/07ef1342-c1a2-11e4-9a04-5bd429305577.png)

We have also tried an analogous model including two "TIA" states (TIA2). This implies a slight better sensitivity (around 2% better).

![TIA2](https://cloud.githubusercontent.com/assets/7307772/6461640/1083f608-c1a2-11e4-8143-f400ab44620d.png)

Finally, we have tested another approach, changing the order of the sates "intron" and "TIA" (using a single "TIA" state). In this case we get a result almost as accurate as with the model TIA2, but using one state less.

![TIA3](https://cloud.githubusercontent.com/assets/7307772/6461641/13f5f5e8-c1a2-11e4-9388-8d0d85e826a6.png)

<a name="objective5"></a>
5. :bar_chart: **Make an assesment of the performance of the model using accuracy measures. Do you find any improvement between models?**

We have already described how to run tests a [couple of sections above](#runtestexample). In order to check the accuracy of each model we calculate the number of true positives (TP), false positives (FP) and false negatives (FN). In this case we can either get the donor position right or wrong, and thus a false positive (we are signaling a position as donor when it isn't) also implies a false negative (we are not signaling the right position as donor). Due to this FP = FN, and the sensitivity **TP/(TP + FN)** is equal to the specificity **TP/(TP + FP)**. A graphic is shown below with the comparison of the sensitivities obtained with each model.

![Accuracies](https://cloud.githubusercontent.com/assets/7307772/6461653/32d7d4ae-c1a2-11e4-8696-b89ab68c483f.png)
