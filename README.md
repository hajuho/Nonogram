# Nonogram

Solver for [Nonogram](https://en.wikipedia.org/wiki/Nonogram)

## Build

```
mkdir build
cd build
cmake ..
make
```

## Run

Run sample
```
./Nonogram
```

Run long sample
```
./Nonogram -l
```

Run sample with showing progress
```
./Nonogram -s
```

Run sample with waiting key input for each step
```
./Nonogram -s -w
```

Run data file
```
./Nonogram ../data/farm_47.txt
```

Run data file with waiting key input for each step
```
./Nonogram -s -w ../data/farm_47.txt
```
