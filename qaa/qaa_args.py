class QAA_ArgumentsAdapter(object):
	def __init__(self, **kwargs):
		for k in kwargs:
			setattr(self, k, kwargs[k])
		pass
	def update(self, **kwargs):
		for k in kwargs:
			# print("Updating: {} <- {}".format(k, kwargs[k]))
			setattr(self, k, kwargs[k])
		pass

