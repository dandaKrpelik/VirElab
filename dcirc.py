import numpy as np
import matplotlib.pyplot as plt
import re
from pyElab import osci



def complex_to_str(a):
	out = 'z is ' +str(a)+'\n\tabs:  ' + str(abs(a)) + '\tdegs:  ' + str(np.rad2deg(np.angle(a)))
	return out
	
def c2str(z):
	out = '{:3g} < {:3g}'.format(abs(z) ,np.rad2deg(np.angle(z)))
	return out
def c2strA(z):
	out = '{:3g}  +j{:.4g}'.format(z.real ,z.imag)
	return out
	
def phaseDiag(xs, labels = None):
	n = len(xs)
	
	if not labels:
		labels = list(range(n))
	
	plt.figure(figsize=(12,10))
	for xi in range(n):
		x = xs[xi]
		plt.plot([0, x.real], [0, x.imag], label = labels[xi], lw = 3)
	plt.legend()
	plt.gca().set_aspect(1)
	plt.plot()
	
class Node:
	_id_counter = 0
	def __init__(self, code='00'):
		self.id = Node._id_counter
		Node._id_counter+= 1
		
		self.typ = code[0]		
		self.tag = code

		if self.typ == 'D': self.open = False

		self.vali = 0
		self.valu = 0
		self.z = 0
		self.neigh = {}
		
		self.loops = []
	
	def __str__(self):
		out = self.tag
		if self.typ == 'D': out += '(+)' if self.open else '(-)'
		out += ': vals= '+str((c2str(self.z), c2str(self.valu), c2str(self.vali)))
		return out

	def str(self, iin):
		out = self.tag
		if self.typ == 'D': out += '(+)' if self.open else '(-)'
		
		u = self.get_u(iin)
		i = self.get_i(iin)
		p = u*np.conjugate(i)
		
		out += ': vals(u,i,p)=  '+c2str(u)+'\t'+c2str(i)+'\t'+c2strA(p)
		return out
	
	def set_val(self, val):
		typ =self.typ
		if typ == 'U': self.valu = val
		elif typ == 'I': self.vali = val
		elif typ == 'D': self.valu = val
		else: self.z = val
		
	def fixedI(self):
		if self.typ == 'I': return True
		if self.typ == 'D':
			if not self.open: return True
		return False		
	
	def fixedU(self):
		if self.typ == 'U': return True
		if self.typ == 'D':
			if self.open: return True
		return False
		
	def get_eq(self, n):
		out = np.zeros(n,dtype = np.complex)
		for loop in self.loops:
			sign = loop.sign(self)
			out[loop.id] += self.z * sign
		return out
		
	def get_i(self, iin):
		i = 0
		for loop in self.loops:
			i += loop.sign(self) * iin[ loop.id ]
		return i

		
	def get_u(self, iin ):
		if self.fixedU():
			return self.valu
		if self.fixedI():
			assert len(self.loops) == 1, self.loops
			loop = self.loops[0]
			
			sgn = 1
			u = 0
			for ni in range(len(loop.nodes)):
				node = loop.nodes[ni]
				if node.id == self.id:
					sgn = loop.sgn[ni]
					u *= sgn
					continue
				u -= sgn * loop.sgn[ni] * node.get_u( iin )
			return u
		
		i = self.get_i(iin)
		u = i*self.z
		return u
		
class Loop:
	def __init__(self, nodes, id = -1):
		self.nodes = nodes
		self.id = id
		self.sgn = [1]*len(nodes)
		
	@property
	def str(self):
		code = ''
		n = len(self.nodes)
		for ni in range(n):
			node = self.nodes[ni]
			sgn = '-' if self.sgn[ni] > 0 else '+' if self.sgn[ni] <0 else 'o'
			code+=sgn+node.tag
		return code
		
		
	def sign(self, node):
		for ni in range(len(self.nodes)):
			if self.nodes[ni].id == node.id:
				return self.sgn[ni]
		return 0
		
	def fixed(self):
		out = None
		for x in self.nodes:
			if x.fixedI(): return True
		return False
		
	def get_fix_i(self):
		for ni in range(len(self.nodes)):
			node = self.nodes[ni]
			if node.fixedI(): return node.vali * self.sgn[ni]
		return 0
		
	def get_eq(self, n):
		outb = 0
		out = np.zeros(n,dtype = np.complex)
		for node in self.nodes:
			if node.fixedU():
				outb -= self.sign(node) * node.valu
				continue
				
			out += self.sign(node) * node.get_eq(n)
		return out, outb
		
class Circuit:
	def __init__(self):
		self.nodes = {}
		self.loops = []
		
		self.sources = {}
		self.diodes = {}
				
	def draw_phase_diag(self, ax, tagsu, tagsi, ism = None):
		if ism is None:
			ism = self.solve_fix()
		
		uus = [ self.nodes[tag].get_u(ism) for tag in tagsu ]
		iis = [ self.nodes[tag].get_i(ism) for tag in tagsi ]
		
		
		maxu = -1
		maxi = -1
		
		for i in range(len(tagsu)):
			node = self.nodes[tagsu[i]]
			x = uus[i]
			maxu = abs(x) if abs(x) > maxu else maxu
			ax.plot( [0,x.real],[0,x.imag] , ':', label = 'U'+tagsu[i])
			
		for i in range(len(tagsi)):
			node = self.nodes[tagsi[i]]
			x = iis[i]
			maxi = abs(x) if abs(x) > maxi else maxi
		
		ICOEF = maxu/maxi/2 if maxu > 0 else 1
		
		for i in range(len(tagsi)):
			node = self.nodes[tagsi[i]]
			x = iis[i]
			if maxu > 0:
				x *= ICOEF
			ax.plot( [0,x.real],[0,x.imag] , label = 'I'+tagsi[i])
			
		ax.scatter([0],[0],color = 'black', s=125, label = '0')
		ax.legend()
		ax.grid('both')
		ax.set_aspect(1)
		
		if ICOEF > 0 and ICOEF != 1:
			plt.title('meritko proudu : {:3.2f} (d) = 1 (A)'.format(ICOEF))
		
		return ICOEF
		
				
	def phase_diag(self, tagsu, tagsi, ism = None):
		
		plt.figure(figsize=(12,10))
		ax = plt.gca()
		self.draw_phase_diag(ax, tagsu, tagsi, ism)
		plt.show(block=False)
			
				
	def set_vals(self, dick):
		for key in dick:
			self.nodes[key].set_val(dick[key])
				
	def add_node(self, node):
		tag = node.tag
		if not tag in self.nodes:
			self.nodes[ tag ] = node
			
			if node.typ in ['U','I']:
				self.sources[tag] = node
			if node.typ in ['D']:
				self.diodes[tag] = node
		
	def solve_bin(self, x):
		assert len(x) == len(self.sources), str(x) + '\n'+ str(sources)
		snames = list(self.sources)
		for si in range(len(snames)):
			sn = snames[si]
			self.sources[sn].set_val( x[si] )		
		
		ds = [self.diodes[x] for x in self.diodes]
		dn = len(ds)
		if dn == 0: return self.solve_fix(x)
		
		for d in ds:
			d.open = False	
		
		stop = False
		while not stop:
			fixi = self.solve_fix(x)
			
			stop = True
			for d in ds:
				u = d.get_u(fixi)
				if d.open and u < d.valu:
					d.open = False
					stop = False
					continue
				if not d.open and u > d.valu:
					d.open = True
					stop = False
					continue
			
		
		u = { x:self.nodes[x].get_u( fixi ) for x in self.nodes }
		
		#print('u',u)
		#print([x.tag+' is '+( 'open' if x.open else 'closed' ) for x in ds])
		return fixi
		
		
	def solve_fix(self, x = None):
		if x:
			assert len(x) == len(self.sources), str(x) + '\n'+ str(sources)
			snames = list(self.sources)
			for si in range(len(snames)):
				sn = snames[si]
				self.sources[sn].set_val( x[si] )
			
		M, b = self.get_M()
		print('M,b', M,b)
		
		n = len(self.loops)
		fixi = np.zeros(n, dtype = np.complex)
		mask = []
		for i in range(n):
			loop = self.loops[i]
			if loop.fixed():
				fixi[i] = loop.get_fix_i()
			else:
				mask+=[i]
				
		b -= M.dot(fixi)
		#print('b', b)
		A = np.matrix(M[mask][:,mask])
		bb = b[mask]
		#print('mask,A,bb', mask,A,bb)
		ism = np.linalg.solve(A,bb)
		#print('ism',ism)
		fixi[mask] = ism
		#print('fixi',fixi)
		
		#u = { x:self.nodes[x].get_u( fixi ) for x in self.nodes }
		#return [ (x,u[x]) for x in self.nodes ]
		return fixi
		
	def tellegen(self, ism):
		for loop in self.loops:
			ss = 0
			for n in loop.nodes:
				ss += loop.sign(n) * n.get_u(ism)
			print('loop ',loop.nodes,' \n\t\t sum P = ',ss)
		
	def get_M(self):
		n = len(self.loops)
		M = np.zeros((n,n), dtype = np.complex)		
		b = np.zeros(n, dtype = np.complex)		
		for i in range(n):
			loop = self.loops[i]
			
			eq = loop.get_eq(n)
			M[i,:] = eq[0]
			b[i] = eq[1]
		
		#print('MB:',M,b)
		return M, b
		
	def add_loop(self, loop = []):
		# loop is list of tags
		# adds edges between nodes in circuit graph
		
		n = len(loop)
				
		nodes = []
		for i in range(n ):
			tag = loop[i][1:]
			sgn = loop[i][0]
			
			if not tag in self.nodes:
				self.add_node( Node( tag ) )
			
			node = self.nodes[ tag ]
			nodes += [node]
		
		loop_obj = Loop( nodes, id = len(self.loops) )
		for i in range(n ):
			tag = loop[i][1:]
			sgn = loop[i][0]
			loop_obj.sgn[i] = 1 if sgn == '-' else -1
			nodes[i].loops += [ loop_obj ]
		
		self.loops += [loop_obj]
		print('added loop ',loop_obj.str)

	

	def get_loop_M(self):
		n = len(self.loops)
		b = np.zeros(n, dtype = np.complex)
		M = np.zeros((n,n), dtype = np.complex)
		
		mask = np.zeros(n,dtype = np.bool)
		
		for i in range(n):
			loopi = self.loops[i]
			for k in range(len(loopi)):
				sgn = 1 if loopi[k][0] == '+' else -1
				tagk = loopi[k][1:]
				node = self.nodes[ tagk ]
				if node.type == 'U':
					b[i] += sgn*node.val
				elif node.type == 'M':
					sgn *= -1
					for j in range(n):
						loopj = self.loops[j]
						jtags = [x[1:] for x in loopj]
						for m in range(len(jtags)):
							if jtags[m] == 'L'+tagk[2:]: #tagk:		## corresponding inductor
								sgnj = sgn if loopj[m][0] == '-' else -sgn
								M[i,j] += sgn * node.val
					
				else:
					sgn *= -1
					if node.type == 'I':
						mask[k] = 1
						continue
					
					for j in range(n):
						loopj = self.loops[j]
						jtags = [x[1:] for x in loopj]
						for m in range(len(jtags)):
							if jtags[m] == tagk:
								sgnj = sgn if loopj[m][0] == '-' else -sgn
								M[i,j] += sgn * node.val
								
		return M, b
		
	def show_phases(self, si, oset = 0.05):
		xi = {key:0 for key in self.nodes}
		xu = {key:0 for key in self.nodes}
		
		for i in range(len(self.loops)):
			loop = self.loops[i]
			for key in loop:
				sgn = 1 if key[0] == '-' else -1
				elem = key[1:]
				
				xi[elem] +=  sgn * si[i]
				
				if not elem[0] == 'U':
					xu[elem] += xi[elem] * self.nodes[elem].val
				else:
					xu[elem] = self.nodes[elem].val
		
		for elem in xi:
			if elem[0] == 'U':
				xi[elem] *= -1
				#xu[elem] *= -1
			
		print(xi)
		print(xu)
		
		plt.figure()
		for elem in xi:
			plt.plot( [0,xi[elem].real] , [0,xi[elem].imag ], label = 'I_'+elem, lw =3, alpha = 0.5, linestyle = ':')
			plt.plot( [0,xu[elem].real] , [0,xu[elem].imag] , label = 'U_'+elem, alpha = 0.5)
		
		for elem in xi:
			plt.text(xi[elem].real, xi[elem].imag+oset, 'I_'+elem)
			plt.text(xu[elem].real, xu[elem].imag-oset, 'U_'+elem)
			
		
		plt.legend()
		plt.gca().set_aspect(1)
		plt.show(block = False)
			
def str2circ( code ):
	print( 'constructing circuit object for:\n\t' +code)
	
	out = Circuit()
	loops = code.split(',')
	for loop in loops:
		print('\t-- loop '+str(loop))
		pos = loop.find('+-')
		
		p = re.compile(r'[+-]\w+')
		tagloop = p.findall( loop )
		out.add_loop(tagloop)
	
	return out
	
	
	
if __name__ == '__main__':
		code = '+U-R1-D1,+U-R1-R2'
		c = str2circ( code )
		
		for name in c.nodes:
			c.nodes[name].set_val(3)
			
