from scipy.optimize import minimize

class FunctionWrapper:
	def __init__(self, f, x, eps, ftol):
		self.f = f
		self.eps = eps
		self.ftol = ftol
		self.start_x = x
		self.last_x = x
		self.grad_started = False
		self.grad_finished = False
		
	def check_diff(self, x1, x2):
		assert len(x1) == len(x2)
		n_diff = 0
		s_diff = 0
		for xx1, xx2 in zip(x1, x2):
			if abs(xx1 - xx2) >= self.PREC:
				n_diff += 1
				s_diff += abs(xx1 - xx2)
		return n_diff, s_diff

	def check_if_step_not_grad(self, x):
		x = list(x)
		n_diff, s_diff = self.check_diff(x, self.last_x)
		if n_diff == 1 and abs(s_diff - eps) <= self.PREC:
			return False
		else:
			self.last_x = x
			return True

	def calculate(self, x):
			t = self.check_if_step_not_grad(x)
			if t:
				if self.grad_started:
					self.grad_finished = True
					return 0
				return self.f(x)
			else:
				self.grad_started = True
				if self.grad_finished:
					return 0
				return self.f(x)
	
	def make_opt_step(self):
		minimize(self.f, self.last_x, args=(), method='SLSQP', jac=None, bounds=[], constraints=[], tol=None, callback=None,
			options={'disp': False, 'eps': self.eps, 'maxiter': 1, 'ftol': ftol})
		return self.last_x
























