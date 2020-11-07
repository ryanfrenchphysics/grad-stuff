import multiprocessing as mp

def func(args):
	first,second = args
	#print(first)
	return first + second

def main(n,n_process):
	pool = mp.Pool(processes = n_process)
	tasks = [(.05*i,.25) for i in range(0,n)]
	results = pool.map(func,tasks)
	print(results)
main(40,8)
