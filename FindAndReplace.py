# Read in the file
with open('StateConstraints.m', 'r') as file :
  filedatain = file.read()

# the dictionary has target_word : replacement_word pairs (mind the order!)
wordDic = {
'vy': 'xi{k}(1)',
'vx': 'xi{k}(2)',
'thetad': 'xi{k}(4)',
'wfl': 'xi{k}(5)',
'wfr': 'xi{k}(6)',
'wrl': 'xi{k}(7)',
'wrr': 'xi{k}(8)',
'lf': 'params.lf',
'lr': 'params.lr',
'rw': 'params.rw',
'w': 'params.w',
'eps10': 'e{k}(10)',
'eps11': 'e{k}(11)',
'eps12': 'e{k}(12)',
'eps13': 'e{k}(13)',
'eps14': 'e{k}(14)',
'eps15': 'e{k}(15)',
'eps16': 'e{k}(16)',
'eps1': 'e{k}(1)',
'eps2': 'e{k}(2)',
'eps3': 'e{k}(3)',
'eps4': 'e{k}(4)',
'eps5': 'e{k}(5)',
'eps6': 'e{k}(6)',
'eps7': 'e{k}(7)',
'eps8': 'e{k}(8)',
'eps9': 'e{k}(9)',
'rparams.w': 'rw',
}

for key in wordDic:
    filedatain = filedatain.replace(key, wordDic[key])

# Write the file out again
with open('StateConstraintsYALMIPformat.m', 'w') as file:
  file.write(filedatain)
