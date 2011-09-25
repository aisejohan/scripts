myN = 20;
myepsilon = 0.0001;

findsmall(f,p,n) = 
{
	local(d,c,i,polygon,j,e,g);
	if(isprime(p),,error("Second argument should be a prime."));
	d = poldegree(f);
	g = x^d;
	polygon = newtonpoly(f,p);
	for(i=0,d-1,
		c = polcoeff(f,i);
		e = floor(sum(j=1,d-i,polygon[d+1-j]));
		c = c/p^e;
		c = lift(Mod(c,p^n));
		if((c>(p^n)/2),c=c-p^n);
		g = g + c*p^e*x^i;
	);
	return(g);
}

check_symmetry(f,p,w) =
{
	local(d,c);
	if(isprime(p),,error("Second argument should be a prime."));
	d = poldegree(f);
	c = polcoeff(f,0)/p^(w*d/2);
	if((c^2 - 1),
		print("Not a Weil polynomial.");
		return(0)
	);
	if((c == 1),
		if(x^d*p^(-w*d/2)*subst(f,x,p^w/x)-f,
			print("Not symmetric.");
			return(0)
		,
			print("The sign is ",(-1)^d,".");
			return(1)
		)
	);
	if((c == -1),
		if(x^d*p^(-w*d/2)*subst(f,x,p^w/x)+f,
			print("Not symmetric.");
			return(0)
		,
			print("The sign is ",-(-1)^d,".");
			return(1)
		)
	);
}

/* Finds Weil polynomial. */
findweil(f,p,initial) =
{
	local(d,polygon,w,n,g,lijst,success);
	if(isprime(p),,error("Second argument should be a prime."));
	if(initial,,initial=1+floor(log(poldegree(f)*p)/log(p)));
	d = poldegree(f);
	polygon = newtonpoly(f,p);
	w = 2*sum(i=1,d,polygon[i])/d;
	print1("The weight is ",w,".");
	n=myN+initial+1;
	while(n>=initial,
		n = n - 1;
		g = findsmall(f,p,n);
		lijst = abs(polroots(g));
		success = 1;
		for(i=1,matsize(lijst)[1],
			if(((p^(w/2)-myepsilon > lijst[i]) ||
				(lijst[i] > p^(w/2)+myepsilon)),
				success = 0
			)
		);
		if(success,
			print1(" Something found for n=",n,". ");
			if(check_symmetry(g,p,w),
				print("Success!");
				return(g)
			)
		)
	);
	print(" No luck this time!");
	return(0);
}

/* Assumes that A times its denominator has positive weight. */
do_it(A)=
{
	local(p,lijst,m,c,l,d,f,g);
	f = charpoly(A);
	d = matdet(denominator(A)*A);
	lijst = factor(d,0);
	l = matsize(lijst)[1];
	m = lijst[1,2];
	c = 1;
	for(i=2,l,if((lijst[i,2] > m),c=i;m=lijst[i,2]));
	p = lijst[c,1];
	print("The prime is ",p,".");
	g = findweil(f,p);
	if(g,,error("Did not work."));
	print("The valuation of f-g is ",valuation(f-g,p));
	print("The valuation of subst(g,x,A) is ",valuation(subst(g,x,A),p));
	print("The factorization of g is ",factor(g));
	return(g);
}

repeat()=
{
	system("./tester | tee uit");
	read("uit");
	g = do_it(B);
	print(g);
	p = factor(abs(subst(g, x, 0)))[1,1];
	newtonpoly(g, p);
}
