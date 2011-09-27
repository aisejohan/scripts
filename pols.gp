myN = 20;
myepsilon = 0.0001;

check_symmetry(f,p,w) =
{
	local(d,c);
	if(isprime(p),,error("Second argument should be a prime."));
	d = poldegree(f);
	c = polcoeff(f,0)/p^(w*d/2);
	if((c^2 - 1),
		return(0)
	);
	if((c == 1),
		if(x^d*p^(-w*d/2)*subst(f,x,p^w/x)-f,
			return(0)
		,
			return(1)
		)
	);
	if((c == -1),
		if(x^d*p^(-w*d/2)*subst(f,x,p^w/x)+f,
			return(0)
		,
			return(1)
		)
	);
}

check_weil(g,p,w) =
{
	local(lijst,success,i);
	lijst = abs(polroots(g));
	for(i=1,matsize(lijst)[1],
		if(((p^(w/2)-myepsilon > lijst[i]) || (lijst[i] > p^(w/2)+myepsilon)), return(0))
	);
        return(1)
}

create_list(p,w,d) =
{
	local(f,g,i,j,a,b,c,N,odd);
	if(w*d % 2,return());
	odd = d % 2;
	l = floor(d/2);
	a = vector(l);
	b = vector(l);
	c = vector(l);
	i = 1;
	N = 1;
	g = x^d + p^(d*(w/2));
	f = x^d - p^(d*(w/2));
	while(i <= l,
		b[i] = floor(binomial(d, i)*p^((w/2)*i));
		N = N*(2*b[i] + 1);
		a[i] = - b[i];
		c[i] = p^((w/2)*(d - 2*i));
		if(i < d/2,
			g = g + a[i]*(x^(d - i) + c[i]*x^i);
			f = f + a[i]*(x^(d - i) - c[i]*x^i)
		,
			/* i = l = d/2 */
			g = g + a[i]*x^i
		);
		i = i + 1
	);
	j = 1;
	while(j <= N,
/*		if(j % 100,,print("Done ",j," out of ", N));
		print(a);
		print(f);
		print(g);
		if(check_symmetry(f,p,w),,error("Not symmetric"));
		if(check_symmetry(g,p,w),,error("Not symmetric")); */
		if(check_weil(g,p,w), print(g));
/*		if(check_weil(f,p,w), print("DUPLICATES   ", f)); */
		if((odd || (b[d/2] == 0)),
			if(check_weil(f,p,w), print(f))
		);
		i = 1;
		while(i <= l,
			if(a[i] + 1 <= b[i],
				if(i < d/2,
					g = g + (x^(d - i) + c[i]*x^i);
					f = f + (x^(d - i) - c[i]*x^i)
				,
					g = g + x^i
				);
				a[i] = a[i] + 1;
				break
			,
				if(i < d/2,
					g = g + (-b[i]-a[i])*(x^(d - i) + c[i]*x^i);
					f = f + (-b[i]-a[i])*(x^(d - i) - c[i]*x^i)
				,
					g = g + (-b[i]-a[i])*x^i
				);
				a[i] = - b[i];
			);
			i = i + 1
		);
		j = j + 1;
	)
}

create_list(2, 2, 6)
quit
