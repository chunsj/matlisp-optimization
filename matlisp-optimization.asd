(asdf:defsystem matlisp-optimization
  :name "matlisp-optimization"
  :version "0.1"
  :author "Akshay Srinivasan <akshays@cs.washington.edu>"
  :licence "LLGPL, see LICENSE"
  :description "Optimization package for Matlisp"
  :depends-on (#:matlisp #:fiveam)
  :serial t
  :components  
  ((:module "package-declaration" :pathname "src" :components ((:file "packages")))
   (:module "src" :depends-on ("package-declaration")
	    :components
	    ((:module "quasi-newton" :components ((:file "bfgs") (:file "sr1")))
	     (:file "linesearch")
	     (:file "bfgs" :depends-on ("quasi-newton" "linesearch"))))))
