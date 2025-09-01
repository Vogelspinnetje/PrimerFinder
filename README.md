# PrimerFinder
## Description
PrimerFinder is een programma die het mogelijk maakt om primers te vinden. Het prgramma zoekt primers in een aangegeven locus in een *genbank* bestand. Het script geeft alle mogelijke primers in aangegeven gebied. Je krijgt zowel forward als backward primers. Als er geen primers gevonden zijn, krijg je een melding. Als de locus niet gevonden wordt, krijg je ook een melding. Het script zorgt er altijd voor dat een primer eindigt met een *Guanine* of *Cytosine* base. 

- De forward primers worden weergeven als 5' - 3' (**geen** reverse complement) ten opzichte van het genbank bestand
- De backward primers worden weergeven als 5' - 3' **reverse complement** ten opzichte van het genbank bestand


## Requirements
Het script is geschreven in `python 3.12.0`. Er wordt aangeraden om minstens `python 3.9` te gerbuiken. Requirements over python packages zijn te vinden in `requirements.txt`.

## Usage
Het is helemaal niet lastig om dit script te gebruiken, je hoeft hem alleen maar te runnen. Het is vooral belangrijk om de instellingen goed in te voeren. Dit doe je in de config.yaml. In deze repo staat een voorbeeld van hoe een config.yaml eruit ziet (met settings die ik zelf heb gebruikt om mijn primers te vinden).