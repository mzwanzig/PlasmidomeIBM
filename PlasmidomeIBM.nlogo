;************************************************************************
;
; PlasmidomeIBM - Simulates a diverse plasmid community
;                 in structured and mixed environments;
;                 considering the modular basis of plasmid
;                 costs and transfer probabilities in relation
;                 to costs associated with the conjugation module
;
; © Martin Zwanzig 2020
;
;   M Zwanzig, U Berger (2020): Modelling the impact of plasmid
;   interactions on coexistence and antibiotic resistance in evolving
;   plasmid communities. eLife (under review; 14.07.2020)
;
; Created under NetLogo Version 6.1.1
;
;************************************************************************

extensions [vid]

globals [
  ;STOP CONDITIONS
  ;stop-when-resistance-is-lost? ;all plasmids confering antibiotic resistance are lost
  ;simulation-time               ;maximum number of ticks to execute

  ;VARIABLES DETERMINING TRAITS AND DIVERSITY OF INITIAL PLASMID POPULATION
  ;initial-plasmid-density
  ;rm-mean
  ;am-max
  ;cm-mean
  ;min-rc
  ;min-cc
  ;ec-mean
  ;dev-strength
  ;inc-numbers

  ;VARIABLES DETERMINING GLOBAL REGIME
  ;seg-prob             ;probability that a plasmid is not transfered to a daughter cell during bacterial fission
  ;mortality            ;probability that a bacterium and all associated plasmids die (=washout/replacement)
  ;mixed-environment?   ;turning complete mixing of bacteria on (= planktonic phase) or off (= surface attached phase)
  ;immigration          ;bacteria may immigrate with a probability of immigration in relation to mortality;
                        ;frequency of plasmid-free and plasmid-bearing bacteria according to definitions for the initial population
  ;baa                  ;bacteriostatic antibiotic action (=growth rate reduction) on non-resistant cells (=those that do not carry the plasmid encoding the resistance)
  ;initial-ARP?         ;defines if and how the initial set of plasmids encoding antibiotic resistancs is created (not at all or with mean or random properties)
  ;conjugative-ARP?     ;defines if the initial ARP plasmids should be conjugative or not
  ;consider-single-plasmid-dynamics? ;the model is initalized only with a single random plasmid from the population
  ;

  ;REPORTER:
  Fc            ;number of plasmid-free bacteria
  Pc            ;number of plasmid-bearing bacteria
  ARB           ;number of resistant bacteria (bearing a plasmid confering resistence to antibiotics)
  plasmid-count ;total number of plasmids
  inc-div       ;number of different plasmid incompatibility groups
  plasmid-div   ;number of different plasmid types
  plasmid-host-div ;number of different plasmid-host-associations
  ARP           ;number of plasmids confering antibiotic resistance
  I1            ;number of plasmids belonging to incompatibility group 1
  I2            ;number of plasmids belonging to incompatibility group 2
  I3            ;number of plasmids belonging to incompatibility group 3
  I4            ;number of plasmids belonging to incompatibility group 4
  I5            ;number of plasmids belonging to incompatibility group 5
  can-pb.       ;plasmid-burden of the candidate plasmid if only single plasmid dynamics is considered
  can-tp.       ;transfer probability of the candidate plasmid if only single plasmid dynamics is considered
  reproduction  ;number of reproduction events per tick
  host-fitness-list ;a list containing the actual sum of plasmid costs for each plasmid-bearing host-bacterium

  ;HELPER VARIABLE
  pid-counter           ;used to give unique plasmid identities to every unique set of plasmid properties
  resistance            ;used to check if any of the plasmids of the focal cell provides resistance to antibiotics
  Fcol                  ;color of plasmid-free bacteria
  Pcol                  ;color of plasmid-bearing bacteria
  vid-counter           ;used to make a frame for the video only every x ticks/generations
  ;evaluate-pb-tp-relations? ;if on, the 'pb-tp-relations' plot is updated every tick, otherwise not (speeds up simulation)
  my-fitness            ;used to calculate host-fitness
]

patches-own [
  op-help. ;operational helper variable
  target.  ;operational helper variable
]

turtles-own [
  pid.     ;plasmid type identity reflecting a unique set of plasmid properties (helps to identify the abundance of this set)
  rm.      ;replication module
  am.      ;accessory module
  cm.      ;conjugation module
  ec.      ;efficiency of conjugation (related to costs of the conjugation module cm.)
  inc.     ;plasmids incompatibility group (a bacterium can harbour only one plasmid of each group at the same time)
  pb.      ;helper variable to store the plasmid-burden resulting from rm. + am. + cm. (reduces the growth of the host bacterium)
  tp.      ;helper variable for plasmids transfer probability (0 < cm. * ec. < 1)
  res.     ;defines if antibiotic resistance genes are present (only for plasmids with am. costs > 0.05)
]

;*****************************RESET**************************************
to reset ; Resets the model parameters to a default setting
  set record-video? false
  set world-dim 250
  set initial-plasmid-density 0.8
  set rm-mean 0.03
  set am-max 0.25
  set cm-mean 0.2
  set min-rc 0.01
  set min-cc 0.01
  set ec-mean 2
  set dev-strength 0.5
  set inc-numbers 3
  set mortality 0.2
  set immigration 0
  set seg-prob 0.001
  set mixed-environment? false
  set consider-single-plasmid-dynamics? false
  set surface-exclusion? true
  set inc-seg-mechanism "random-daughter-load"
  set baa 0
  set initial-ARP? "none"
  set conjugative-ARP? false
  set ARP-prop 0
  set generations_of_antibiotic_presence 0
  set simulation-time 10000
  set stop-when-resistance-is-lost? false
  set evaluate-tp-pb-relations? true
end

;*****************************SETUP**************************************

to setup
  clear-all
  vid:reset-recorder

  resize-world 0 (world-dim - 1) 0 (world-dim - 1)
  set-patch-size 2 * (250 / world-dim)

  set Fcol white
  set Pcol 9

  ask patches [
    ifelse random-float 1.0 < (1 - mortality)
    [
      set pcolor Fcol
    ]
    [ set pcolor black ]
  ]

  let plasmid-number count patches with [pcolor = Fcol] * initial-plasmid-density
  ask n-of plasmid-number patches with [pcolor = Fcol] [
    create-a-new-plasmid
  ]

  if any? turtles [
    ifelse consider-single-plasmid-dynamics?
    [
      let candidate one-of turtles
      set can-pb. [pb.] of candidate
      set can-tp. [tp.] of candidate
      ask candidate [
        if initial-ARP? != "none" [set res. 1 set color orange]
        let dismissed other turtles
        ask dismissed [
          let dismissed-cell patch-here
          ask candidate [hatch 1[move-to dismissed-cell]]
          die
        ]
      ]
    ]
    [
      if initial-ARP? = "ARP-with-mean-properties" [specify-ARP]
      if initial-ARP? = "ARP-with-random-properties" [select-random-ARP]
      if initial-ARP? = "many-random-ARP-with-nearly-mean-am." [select-many-random-ARP]
    ]
  ]

  reset-ticks
  report-system-state
  do-plots
  draw-tp-pb-relations

  if record-video? [ vid:start-recorder vid:record-interface ]
  make-video
end

to create-a-new-plasmid
  sprout 1 [
    set pid. pid-counter + 1
    set pid-counter pid-counter + 1
    set rm. random-normal rm-mean (rm-mean * dev-strength)
    while [rm. > 1 or rm. < min-rc] [ set rm. random-normal rm-mean (rm-mean * dev-strength) ]
    set am. random-float am-max
    set res. 0
    ifelse random 2 = 0
    [ ; create non-transmissible plasmid
      set cm. 0
    ]
    [ ; create conjugative plasmid
      set cm. random-normal cm-mean (cm-mean * dev-strength)
      while [cm. > 1 or cm. < min-cc] [ set cm. random-normal cm-mean (cm-mean * dev-strength) ]
    ]
    set pb. rm. + am. + cm.
    set ec. random-normal ec-mean dev-strength
    while [ec. < 0] [ set ec. random-normal ec-mean dev-strength ]
    set tp. cm. * ec.
    while [tp. > 1 or tp. < 0] [
      set ec. random-normal ec-mean dev-strength
      while [ec. < 0] [ set ec. random-normal ec-mean dev-strength ]
      set tp. cm. * ec.
    ]
    set inc. 1 + random inc-numbers
    set shape "line"
    if inc. = 1 [ set color cyan ]
    if inc. = 2 [ set color sky ]
    if inc. = 3 [ set color blue ]
    if inc. = 4 [ set color violet ]
    if inc. = 5 [ set color magenta ]
    if inc. > 5 [ set color 15 + 10 * inc. ]
    ;set color 15 + (inc. * 10)
    set pcolor Pcol
    set heading 0 + pb. * 90
  ]
end

to specify-ARP
  let center round (world-dim / 2)
  ask patch center center [ ask one-of turtles-on neighbors
    [
      set res. 1
      set color orange
      set rm. rm-mean
      set am. am-max / 2
      ifelse conjugative-ARP? [set cm. cm-mean] [set cm. 0]
      set pb. rm. + am. + cm.
      set ec. ec-mean
      set tp. cm. * ec.
      set inc. 1
      let targets. patches in-radius 5 with [pcolor != black]
      ask other turtles-on targets. [die]
      hatch (count targets. - 1) [
        move-to one-of targets.
        while [any? other turtles-here]
        [move-to one-of targets.]
      ]
    ]
  ]
end

to select-random-ARP
  let center round (world-dim / 2)
  ask patch center center [
    ifelse conjugative-ARP?
    [
      ask one-of turtles in-radius 3 with [cm. > 0]
      [
        set res. 1
        set color orange
        set inc. 1
        let targets. patches in-radius 5 with [pcolor != black]
        ask other turtles-on targets. [die]
        hatch (count targets. - 1) [
          move-to one-of targets.
          while [any? other turtles-here]
          [move-to one-of targets.]
        ]
      ]
    ]
    [
      ask one-of turtles in-radius 3 with [cm. = 0]
      [
        set res. 1
        set color orange
        set inc. 1
        let targets. patches in-radius 5 with [pcolor != black]
        ask other turtles-on targets. [die]
        hatch (count targets. - 1) [
          move-to one-of targets.
          while [any? other turtles-here]
          [move-to one-of targets.]
        ]
      ]
    ]
  ]
end

to select-many-random-ARP
  let mean-am mean [am.] of turtles
  let sd-am standard-deviation [am.] of turtles
  ask n-of round (ARP-prop * count turtles) turtles with [(mean-am - sd-am) < am. and am. < (mean-am + sd-am)]
  [ set res. 1 set color orange - random-float 0.5 + 0.25 ]
end

;*******************************GO***************************************

to go
  if simulation-time != 0 and ticks >= simulation-time [print "tick max. reached" stop]
  if (initial-plasmid-density > 0 and Pc = 0) [print "plasmid extinction" stop]
  if stop-when-resistance-is-lost? [ if (ARP = 0) [print "resistance lost" stop] ]
  set reproduction 0 ;reset to track the number of reproduction events per time step
  if generations_of_antibiotic_presence > 0 [ if ticks >= generations_of_antibiotic_presence [set baa 0]]
  repeat (world-width * world-height) * (1 / mortality) ;means repetitions for 1 generation of bacteria
  [
    ;update-asynchronously
    ask one-of patches [
      if (pcolor != black) [

        let resource-availability (count neighbors with [pcolor = black] / 8)
        let total-plasmid-burden sum([pb.] of turtles-here)
        set resistance 0
        if (sum([res.] of turtles-here) > 0) [set resistance 1]
        let fission-probability ((1 - total-plasmid-burden) * resource-availability * (1 - baa * (1 - resistance)))
        let total-transfer-probability sum([tp.] of turtles-here)
        let transfer-probability (total-transfer-probability * resource-availability)

        let random-number random-float 1

        ifelse (random-number < mortality) [ lyse ]
        [
          ifelse (random-number < mortality + fission-probability) [ perform-fission ]
          [
            if (random-number < mortality + fission-probability + transfer-probability) [ attempt-transfer ]
          ]
        ]

        if immigration > 0 ;immigration causes that a lattice cell is replaced by an immigrant
        [
          if (random-float 1 < immigration) [ immigrate ]
        ]
      ]
    ]
  ]
  tick
  report-system-state
  do-plots
  if evaluate-tp-pb-relations? [ draw-tp-pb-relations ]
  ;export-gui
  make-video
end

to make-video
  if record-video? [
    if vid-counter = 0
    [ vid:record-interface
      set vid-counter 5 ] ; CHANGE HERE, e.g. to 0 in order to make a frame for the video every tick
    set vid-counter vid-counter - 1
  ]
end

to save-recordings
  if record-video? [ vid:save-recording "out.mp4" ]
end

to lyse ;patch-procedure
  set pcolor black
  ask turtles-here [ die ]
end

to perform-fission ;patch-procedure
  ifelse mixed-environment?
  [ set target. one-of patches with [pcolor = black] ]
  [ set target. one-of neighbors with [pcolor = black] ]
  if pcolor = Pcol [
    ; here two different mechanisms can be used or tested; there are only differences when the mother cell carries more than
    ; one plasmid type of the same incompatibility group (when entry-exclusion is off!!!); In this case one of the incompatible plasmid
    ; types is eliminated either in both cells ("identical-daughter-load" mechanism) or randomly in each cell ("random-daughter-load" mechanism),
    ; so that the arising daughter cells can carry different plasmid types of the same incomaptibility group (or by chance both
    ; have the same and one type is completely eliminated from both cells as in the "identical-daughter-load" mechanism);
    ; in reality this process depends on the plasmid copy number and some other mechanisms preventing missegregation that we do not
    ; explicitly consider, but that could prolong the presence of incompatible plasmid in the same bacterial cell
    ; (our model does not prevent the co-occurrence of incompatible plasmids in the same bacteria cell, but allows this only for one generation)
    if inc-seg-mechanism = "random-daughter-load"
    [ ; in this version each daughter cell has a random plasmid load originating from the mother cells plasmid load
      ; segregation can additionally remove a plasmid from one of both daughter cells
      let source. self
      ask turtles-here [
        hatch 1
        [
          move-to target.
          if random-float 1 <= seg-prob [
            ask one-of (turtle-set self myself) [die]
          ]
        ]
      ]
      ask turtles-here [ if any? other turtles-here with [inc. = [inc.] of myself] [die] ] ;randomly remove plasmids with the same inc-group (only one remains in each daughter cell)
      ask turtles-on target. [ if any? other turtles-here with [inc. = [inc.] of myself] [die] ] ;randomly remove plasmids with the same inc-group (only one remains in each daughter cell)
    ]
    if inc-seg-mechanism = "identical-daughter-load"
    [ ; in this version daughter cells have identical plasmid loads (random incompatible plasmids are removed from both cells)
      ; only segregation can lead to a difference, since this can remove a plasmid only from one of both daughter cells
      ask turtles-here [
        let target-plasmids turtles-on target.
        ifelse not any? target-plasmids with [inc. = [inc.] of myself]
        [
          if random-float 1 > seg-prob [hatch 1 [ move-to target. ]]
        ]
        [ die ]
      ]
    ]
  ]
  ask target. [ set pcolor [pcolor] of myself ]
  if count turtles-here = 0 [set pcolor Fcol]
  if count turtles-on target. = 0 [ask target. [set pcolor Fcol]]
  set reproduction reproduction + 1
end

to attempt-transfer ;turtle-procedure
  ifelse mixed-environment?
  [
    set target. one-of patches
  ]
  [
    set target. one-of neighbors
  ]
  if [pcolor != black] of target. [
    let total-tp sum ([tp.] of turtles-here)
    let random-number random-float total-tp
    set op-help. 0
    ask turtles-here [
      if (random-number < tp. + op-help.) [ ; processes that only one of the plasmids can be transfered
        let target-plasmids turtles-on target.
        if not any? target-plasmids with [pid. = [pid.] of myself][ ;prevents the entry of another, identical plasmid type into the recipient
          ifelse surface-exclusion?
          [ ; prohibits transfer, when the recipient already carries a plasmid of the same incompatibility group
            ; (incompatibility after cell fission is obsolot, since no cell can carry more than one plasmid type of the same inc-group)
            if not any? target-plasmids with [inc. = [inc.] of myself]
            [
              hatch 1 [
                move-to target.
              ]
              ask target. [ set pcolor [pcolor] of myself ]
            ]
          ]
          [ ; transfer happens no matter what incompatibility groups are present
            ; (incompatibility then takes effect after cell fission)
            hatch 1 [
              move-to target.
            ]
            ask target. [ set pcolor [pcolor] of myself ]
          ]
        ]
      ]
      set op-help. op-help. + tp.
    ]
  ]
end

to immigrate
  ask one-of patches [
    ask turtles-here [die]
    ifelse initial-plasmid-density = 0
    [
      ifelse random-float 1 < 0.5
      [ create-a-new-plasmid ]
      [ set pcolor Fcol ]
    ]
    [
      ifelse random-float 1 < initial-plasmid-density
      [ create-a-new-plasmid ]
      [ set pcolor Fcol ]
    ]
  ]
end

;****************************PLOTTING************************************

to report-system-state
  set Fc count patches with [pcolor = Fcol]
  set Pc count patches with [pcolor = Pcol]

  set plasmid-count count turtles

  if plasmid-count > 1 [
    set inc-div length remove-duplicates [inc.] of turtles
    set plasmid-div length remove-duplicates [pid.] of turtles
    report-plasmid-host-diversity
    report-inc-number-frequencies
    report-resistances
    report-host-fitness
  ]
end

to report-plasmid-host-diversity
  ask patches [set op-help. sum [pb.] of turtles-here]
  set plasmid-host-div length remove-duplicates [op-help.] of patches with [op-help. > 0] ;means (sum of pb > 0) -> only plasmid-bearers are considered
  ask patches [set op-help. 0]
end

to report-inc-number-frequencies
  if inc-numbers > 1 [
    set I1 count turtles with [inc. = 1]
    set I2 count turtles with [inc. = 2]
    if inc-numbers > 2 [
      set I3 count turtles with [inc. = 3]
      if inc-numbers > 3 [
        set I4 count turtles with [inc. = 4]
        if inc-numbers > 4 [
          set I5 count turtles with [inc. = 5]
        ]
      ]
    ]
  ]
end

to report-resistances
  ask patches [set op-help. 0]
  let ARPs turtles with [res. = 1]
  set ARP count ARPs
  ask ARPs [set op-help. 1]
  set ARB count patches with [op-help. = 1]
  ask patches [set op-help. 0]
end

to report-host-fitness
  set host-fitness-list (list)
  ask patches with [pcolor = Pcol] [
    set my-fitness (1 - sum [pb.] of turtles-here) ; calculates the sum of the costs of all its residing plasmids
    if my-fitness < 0 [ set my-fitness 0 ]
    set host-fitness-list lput my-fitness host-fitness-list
  ]
end

to do-plots
  set-current-plot "plasmid-frequency"
    set-current-plot-pen "total"
    plot plasmid-count
  set-current-plot-pen "ARP"
    plot ARP
    if inc-numbers > 1 [
      set-current-plot-pen "Inc1"
      plot count turtles with [inc. = 1]
      set-current-plot-pen "Inc2"
      plot I2
      if inc-numbers > 2 [
        set-current-plot-pen "Inc3"
        plot I3
        if inc-numbers > 3 [
          set-current-plot-pen "Inc4"
          plot I4
          if inc-numbers > 4 [
            set-current-plot-pen "Inc5"
            plot I5
  ]]]]

  set-current-plot "trait-proportion"
  set-current-plot-pen "F"
  plot Fc / (Fc + Pc)
  set-current-plot-pen "P"
  plot Pc / (Fc + Pc)
  set-current-plot-pen "ARB"
  plot ARB / (Fc + Pc)

  set-current-plot "diversity"
  let div-max max list plasmid-div plasmid-host-div
  if div-max != 0 and plot-y-max > 20 * div-max [set-plot-y-range 0 10 * div-max]
    set-current-plot-pen "plasmids"
    plot plasmid-div
    set-current-plot-pen "plasmid-host-pairs"
    plot plasmid-host-div

  set-current-plot "inc. diversity"
    set-current-plot-pen "inc-div"
    plot inc-div

  set-current-plot "plasmid-burden-frequency"
    set-current-plot-pen "pb"
    set-plot-pen-interval 0.01
    histogram [pb.] of turtles

  set-current-plot "rm-distribution"
    set-current-plot-pen "rm"
    set-plot-x-range 0 rm-mean * 3
    set-plot-pen-interval 0.01
    histogram [rm.] of turtles

  set-current-plot "am-distribution"
    set-current-plot-pen "am"
    set-plot-x-range 0 am-max
    set-plot-pen-interval 0.01
    histogram [am.] of turtles

  set-current-plot "cm-distribution"
    set-current-plot-pen "cm"
    set-plot-x-range 0 cm-mean * 3
    set-plot-pen-interval 0.01
    histogram [cm.] of turtles

  set-current-plot "ec-distribution"
    set-current-plot-pen "ec"
    set-plot-x-range 0 ec-mean * 3
    set-plot-pen-interval 0.01
    histogram [ec.] of turtles

  set-current-plot "Host-fitness"
    set-current-plot-pen "host-fitness"
    set-plot-x-range 0 1
    let min-hfl min host-fitness-list
    if (min-hfl < 0) [ set-plot-x-range min-hfl 1 ]
    set-plot-pen-interval 0.01
    histogram host-fitness-list

end

to draw-tp-pb-relations
  set-current-plot "tp-pb-relation"
  clear-plot
  set-current-plot-pen "I1"
  ask turtles with [inc. = 1] [plotxy [tp.] of self (1 - [pb.] of self)]
  set-current-plot-pen "I2"
  ask turtles with [inc. = 2] [plotxy [tp.] of self (1 - [pb.] of self)]
  set-current-plot-pen "I3"
  ask turtles with [inc. = 3] [plotxy [tp.] of self (1 - [pb.] of self)]
  set-current-plot-pen "I4"
  ask turtles with [inc. = 4] [plotxy [tp.] of self (1 - [pb.] of self)]
  set-current-plot-pen "I5"
  ask turtles with [inc. = 5] [plotxy [tp.] of self (1 - [pb.] of self)]
  set-current-plot-pen "other"
  ask turtles with [inc. > 5] [plotxy [tp.] of self (1 - [pb.] of self)]
  set-current-plot-pen "ARP"
  ask turtles with [res. = 1] [plotxy [tp.] of self (1 - [pb.] of self)]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Export-Individual-Level-Data;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-file
  file-close-all
  if file-exists? (word behaviorspace-experiment-name behaviorspace-run-number "_ind-level.txt")
  [ file-delete (word behaviorspace-experiment-name behaviorspace-run-number "_ind-level.txt") ]
  file-open (word behaviorspace-experiment-name behaviorspace-run-number "_ind-level.txt")
  file-write "behaviorspace-run-number"
  file-write "ticks"
  file-write "pid"
  file-write "pid.count"
  file-write "xcor"
  file-write "ycor"
  file-write "pb."
  file-write "rm."
  file-write "cm."
  file-write "am."
  file-write "ec."
  file-write "tp."
  file-write "inc."
  file-write "res."
  file-print ""
end

to check-individual-level-data-export
  if ticks = 0 or ticks = 20 or ticks = 200 or ticks = 2000 [
    export-individual-level-data
  ]
end

to export-gui
  if ticks = 0 or ticks = 10 or ticks = 100 or ticks = 1000 or ticks = 10000 or ticks = 10010 or ticks = 10100 or ticks = 11000 or ticks = 20000 [
    export-view (word "view" ticks ".png")
    export-interface (word "interface" ticks ".png")
    ;export-individual-level-data
  ]
end

to export-individual-level-data
    let unique-pid remove-duplicates [pid.] of turtles
    foreach unique-pid [ x ->
      let spec-pids turtles with [pid. = x]
      let pid-count count spec-pids
      ask one-of spec-pids
      [
        file-write behaviorspace-run-number
        file-write ticks
        file-write pid.
        file-write pid-count
        file-write xcor
        file-write ycor
        file-write pb.
        file-write rm.
        file-write cm.
        file-write am.
        file-write ec.
        file-write tp.
        file-write inc.
        file-write res.
        file-print ""
      ]
    ]
end
@#$#@#$#@
GRAPHICS-WINDOW
306
13
814
522
-1
-1
2.0
1
10
1
1
1
0
1
1
1
0
249
0
249
1
1
1
generations
30.0

BUTTON
13
10
68
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
71
10
126
43
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
13
662
181
695
mortality
mortality
0.01
1
0.2
0.01
1
[-]
HORIZONTAL

SLIDER
6
219
176
252
initial-plasmid-density
initial-plasmid-density
0
1
0.8
0.01
1
[-]
HORIZONTAL

SLIDER
8
526
133
559
inc-numbers
inc-numbers
1
20
3.0
1
1
NIL
HORIZONTAL

PLOT
831
138
1229
345
plasmid-frequency
generations
frequency
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"total" 1.0 0 -6459832 true "" ""
"ARP" 1.0 0 -955883 true "" ""
"Inc1" 1.0 0 -11221820 true "" ""
"Inc2" 1.0 0 -13791810 true "" ""
"Inc3" 1.0 0 -13345367 true "" ""
"Inc4" 1.0 0 -8630108 true "" ""
"Inc5" 1.0 0 -5825686 true "" ""

PLOT
1240
503
1439
623
plasmid-burden-frequency
plasmid-burden
frequency
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"pb" 1.0 1 -16777216 true "" ""

TEXTBOX
148
528
301
570
number of incompatibility groups making up the initial plasmid population
11
0.0
1

TEXTBOX
192
657
306
675
washout probability
11
0.0
1

TEXTBOX
185
222
276
250
initial plasmid density
11
0.0
1

SWITCH
320
560
487
593
mixed-environment?
mixed-environment?
1
1
-1000

TEXTBOX
321
527
490
555
simulate 'biofilm phase' (off) or 'planktonic phase' (on)
11
0.0
1

TEXTBOX
183
754
322
822
bacteria may immigrate with a probability of immigration. Plasmid properties are drawn as for the initial pool.
11
0.0
1

TEXTBOX
7
190
271
218
VARIABLES DETERMINING TRAITS AND DIVERSITY OF INITIAL PLASMID POPULATION
11
0.0
1

TEXTBOX
13
575
252
593
VARIABLES DETERMINING GLOBAL REGIME
11
0.0
1

PLOT
1070
531
1230
687
tp-pb-relation
tp
pb
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"other" 1.0 2 -16777216 true "" ""
"I1" 1.0 2 -11221820 true "" ""
"I2" 1.0 2 -13791810 true "" ""
"I3" 1.0 2 -13345367 true "" ""
"I4" 1.0 2 -8630108 true "" ""
"I5" 1.0 2 -5825686 true "" ""
"ARP" 1.0 2 -955883 true "" ""

SLIDER
7
487
131
520
dev-strength
dev-strength
0
1
0.5
0.01
1
NIL
HORIZONTAL

TEXTBOX
148
492
286
517
coefficient of variation for normal distributions
11
0.0
1

SLIDER
6
256
135
289
rm-mean
rm-mean
0.01
0.09
0.03
0.01
1
NIL
HORIZONTAL

SLIDER
7
292
135
325
am-max
am-max
0
1
0.25
0.01
1
NIL
HORIZONTAL

SLIDER
7
328
136
361
cm-mean
cm-mean
0
1
0.2
0.01
1
NIL
HORIZONTAL

TEXTBOX
146
256
296
365
costs for replication module (rm), adaptive module (am) and conjugation module (cm) are initially drawn from normal (mean for rm, cm) or uniform (max for am) distributions.
11
0.0
1

SLIDER
5
450
132
483
ec-mean
ec-mean
1
10
2.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
147
442
305
488
initial mean efficiency of conjugation (relation of cost to transfer probability)
11
0.0
1

PLOT
1239
133
1439
253
rm-distribution
rm.
frequency
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"rm" 1.0 1 -16777216 true "" ""

PLOT
1240
255
1440
375
am-distribution
am.
frequency
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"am" 1.0 1 -16777216 true "" ""

PLOT
1241
377
1439
497
cm-distribution
cm.
frequency
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"cm" 1.0 1 -16777216 true "" ""

PLOT
1239
628
1439
748
ec-distribution
ec.
frequency
0.0
3.0
0.0
10.0
true
false
"" ""
PENS
"ec" 1.0 1 -16777216 true "" ""

SLIDER
7
368
128
401
min-rc
min-rc
0
1
0.01
0.01
1
NIL
HORIZONTAL

TEXTBOX
146
374
296
404
minimal replication module costs
11
0.0
1

SLIDER
7
404
128
437
min-cc
min-cc
0
1
0.01
0.01
1
NIL
HORIZONTAL

TEXTBOX
147
407
297
437
minimal conjugation module costs
11
0.0
1

SLIDER
524
558
696
591
baa
baa
0
1
0.0
0.01
1
[-]
HORIZONTAL

TEXTBOX
705
546
841
601
bacteriostatic antibiotic action - factor for growth rate reduction of sensitive cells
11
0.0
1

PLOT
834
10
1230
132
trait-proportion
generations
frequency
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"F" 1.0 0 -7500403 true "" ""
"P" 1.0 0 -6459832 true "" ""
"ARB" 1.0 0 -955883 true "" ""

SWITCH
524
650
703
683
conjugative-ARP?
conjugative-ARP?
1
1
-1000

INPUTBOX
12
698
173
758
seg-prob
0.001
1
0
Number

TEXTBOX
185
689
331
746
probability that a plasmid is not transfered to one of the daughter cells during bacterial fission
11
0.0
1

TEXTBOX
707
610
844
736
the model can be initialized with a small cluster of bacteria carrying a plasmid-type that confers resistance or with a proportion of random, diverse plasmid types conferring resistence
11
0.0
1

TEXTBOX
859
663
966
681
STOP CONDITIONS
11
0.0
1

INPUTBOX
14
762
175
822
immigration
0.0
1
0
Number

PLOT
834
398
1225
528
diversity
generations
count
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"plasmids" 1.0 0 -13840069 true "" ""
"plasmid-host-pairs" 1.0 0 -2064490 true "" ""

MONITOR
833
349
1019
394
# diverse plasmid types
plasmid-div
0
1
11

MONITOR
1026
349
1230
394
# diverse plasmid-host-pairs
plasmid-host-div
0
1
11

CHOOSER
524
602
702
647
initial-ARP?
initial-ARP?
"none" "ARP-with-mean-properties" "ARP-with-random-properties" "many-random-ARP-with-nearly-mean-am."
0

INPUTBOX
858
679
1019
739
simulation-time
10000.0
1
0
Number

SWITCH
1026
750
1263
783
stop-when-resistance-is-lost?
stop-when-resistance-is-lost?
1
1
-1000

TEXTBOX
942
740
1027
762
0 = unlimited
11
0.0
1

PLOT
857
532
1040
652
inc. diversity
generations
#
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"inc-div" 1.0 0 -16777216 true "" ""

INPUTBOX
13
597
73
657
world-dim
250.0
1
0
Number

CHOOSER
343
706
517
751
inc-seg-mechanism
inc-seg-mechanism
"random-daughter-load" "identical-daughter-load"
0

SWITCH
94
51
228
84
record-video?
record-video?
1
1
-1000

BUTTON
81
122
233
155
NIL
save-recordings
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
523
688
699
721
ARP-prop
ARP-prop
0
1
0.0
0.001
1
NIL
HORIZONTAL

SWITCH
320
600
513
633
consider-single-plasmid-dynamics?
consider-single-plasmid-dynamics?
1
1
-1000

TEXTBOX
1274
754
1433
785
the model stops anyway when all plasmids are lost
11
0.0
1

INPUTBOX
704
742
890
802
generations_of_antibiotic_presence
0.0
1
0
Number

SWITCH
349
666
515
699
surface-exclusion?
surface-exclusion?
0
1
-1000

TEXTBOX
330
755
696
829
incompatibility either takes place by surface exclusion (preventing the entry of incompatible plasmids to the recipient already harboring a plasmid type of the same group) or during cell fission by one of the inc-seg-mechanisms (that remove incompatible plasmids from the daughter cells)
11
0.0
1

TEXTBOX
396
635
510
662
creates only a single plasmid type
11
0.0
1

TEXTBOX
15
68
80
86
1. turn on:
11
0.0
1

TEXTBOX
14
91
238
118
2. press 'setup'-button (top left)\n3. press 'go'-button and run simulation
11
0.0
1

TEXTBOX
14
123
69
141
4. press:
11
0.0
1

TEXTBOX
29
160
243
178
saved as \"out.mp4\" in local directory
11
0.0
1

TEXTBOX
524
527
816
545
Specifications for antibiotics and resistances:
12
0.0
1

TEXTBOX
899
762
1012
818
-> compare to total simualtion-time (consider regime-shift?)
11
0.0
1

SWITCH
1024
692
1234
725
evaluate-tp-pb-relations?
evaluate-tp-pb-relations?
0
1
-1000

MONITOR
1106
465
1198
510
generations
ticks
0
1
11

PLOT
1236
10
1440
130
Host-fitness
remaining growth potential
frequency
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"host-fitness" 1.0 1 -16777216 true "" ""

TEXTBOX
15
52
93
70
Movie maker:
11
0.0
1

TEXTBOX
132
10
231
41
1. press 'setup'\n2. press 'go'
12
0.0
1

TEXTBOX
88
597
257
653
edge length of the quadratic grid as a specification of the model world size (automatic visual scaling applies)
11
0.0
1

BUTTON
237
10
292
43
NIL
reset
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
230
45
306
82
press to reset\nto default para-\nmeterization
9
0.0
1

@#$#@#$#@
## HOW TO RUN A SIMULATION

This is a NetLogo model that can be easily executed using the buttons in the user interface.

Pressing 'setup' initializes the model world, e.g. bacteria and plasmids are generated and distributed over the model world according to the rules determined as a procedure 'setup' in the model code and according to the specifications given by the slider and switches in the user interface. If, for example, antibiotic resistance is to be considered, settings must be made for the number ('ARP-prop') and characteristics of resistance-carrying plasmids ('initial-ARP?', 'conjugative-ARP?'). By default, the model is initialized without antibiotic resistance plasmids.

Pressing 'go' starts a simulation, i.e. all rules that are determined in the procedure named 'go' in the model code are repeatedly applied until the model is forced to stop, e.g. because the specified simulation time is reached or 'go' is pressed again. A number of global variables determine the environmental conditions and should be determined before "setup" is pressed, but they can also be changed for a running model to observe the associated effects, e.g. to study how an increasing or decreasing bacteriostatic antibiotic action ('baa') affects population dynamics.

## CREDITS AND REFERENCES

M Zwanzig, U Berger (2020), Modelling the impact of plasmid interactions on coexistence and antibiotic resistance in evolving plasmid communities. eLife (under review; 14.07.2020)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="default-experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="dynamics-experiment" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>export-individual-level-data
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ARP-experiment" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>check-individual-level-data-export
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;many-random-ARP-with-nearly-mean-am.&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="dynamics-experiment2" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>export-individual-level-data
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="dynamics-experiment3" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>export-individual-level-data
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="dynamics-experiment4" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>export-individual-level-data
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="dynamics-experiment5" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>export-individual-level-data
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ARP-experiment2" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>check-individual-level-data-export
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;many-random-ARP-with-nearly-mean-am.&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ARP-experiment3" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>check-individual-level-data-export
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;many-random-ARP-with-nearly-mean-am.&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ARP-experiment4" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>check-individual-level-data-export
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;many-random-ARP-with-nearly-mean-am.&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ARP-experiment5" repetitions="20" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>check-individual-level-data-export
go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;many-random-ARP-with-nearly-mean-am.&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="20000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="no-surface-exclusion-experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup
setup-file</setup>
    <go>go</go>
    <final>export-individual-level-data
file-close-all</final>
    <metric>Fc</metric>
    <metric>Pc</metric>
    <metric>plasmid-count</metric>
    <metric>inc-div</metric>
    <metric>plasmid-div</metric>
    <metric>plasmid-host-div</metric>
    <enumeratedValueSet variable="world-dim">
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-plasmid-density">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rm-mean">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="am-max">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cm-mean">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-rc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-cc">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ec-mean">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-strength">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-numbers">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seg-prob">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immigration">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mixed-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inc-seg-mechanism">
      <value value="&quot;random-daughter-load&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surface-exclusion?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consider-single-plasmid-dynamics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="baa">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ARP?">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conjugative-ARP?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ARP-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="generations_of_antibiotic_presence">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-time">
      <value value="30000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-when-resistance-is-lost?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-video?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
