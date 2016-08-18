import numpy as np

class ElectricalField(object):
   def __init__(self, electrical_field, X, Y, distance):
      self._electrical_field = electrical_field.copy()
      self._X = X.copy()
      self._Y = Y.copy()
      self._distance = distance

   def X(self):
      return self._X.copy()

   def Y(self):
      return self._Y.copy()

   def distance(self):
      return self._distance

   def electrical_field(self):
      return self._electrical_field.copy()

   def intensity(self):
      return np.sum(np.abs(self.electrical_field())**2,axis=2)
