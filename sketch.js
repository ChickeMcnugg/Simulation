class Fluid {
  constructor(n, scaled, dt, diff, visc, iter) {
    // n - width of simulation
    // dt - time step for calculation
    // diff - diffusion rate
    // visc - viscosity
    // iter - number of iterations for Gauss-Seidel method
    // scaled - scaling factor for drawing to canvas
    
    this.n = n;
    this.dt = dt;
    this.diff = diff;
    this.visc = visc;
    this.iter = iter;
    this.scaled = scaled;
    
    // s - future density
    // density - current density
    // Vx - current X velocities
    // Vy - current Y velocities
    // Vx0 - previous X velocities
    // Vy0 - previous Y velocities
    
    this.s = new Array(n * n).fill(0);
    this.density = new Array(n * n).fill(0);
    this.Vx = new Array(n * n).fill(0);
    this.Vy = new Array(n * n).fill(0);
    this.Vx0 = new Array(n * n).fill(0);
    this.Vy0 = new Array(n * n).fill(0);
  }
  
  step() {
    this.diffuse(1, this.Vx0, this.Vx, this.visc);
    this.diffuse(2, this.Vy0, this.Vy, this.visc);
    
    this.project(this.Vx0, this.Vy0, this.Vx, this.Vy);
    
    this.advect(1, this.Vx, this.Vx0, this.Vx0, this.Vy0);
    this.advect(2, this.Vy, this.Vy0, this.Vx0, this.Vy0);
    
    this.project(this.Vx, this.Vy, this.Vx0, this.Vy0);
    
    this.diffuse(0, this.s, this.density, this.diff);
    this.advect(0, this.density, this.s, this.Vx, this.Vy);
  }
  
  renderD() {
    colorMode(HSB);
    
    for (let i = 0; i < this.n; i++) {
      for (let j = 0; j < this.n; j++) {
        //Fade density
        let d = this.density[this.getIndex(i, j)];
        d = constrain(d - 0.12, 0, 255);
        this.density[this.getIndex(i, j)] = d;
        
        //Draw density map
        fill(255 - d, 255, 255);
        noStroke();
        square(i * scaled, j * scaled, scaled);
      }
    }
  }
  
  renderV() {
    colorMode(HSB);
    
    for (let i = 0; i < this.n; i++) {
      for (let j = 0; j < this.n; j++) {
        //Draw velocity map
        fill(255 - (30 * pow(pow(this.Vx[this.getIndex(i, j)], 2) + pow(this.Vy[this.getIndex(i, j)], 2), 0.5)), 255, 255);
        noStroke();
        square(i * scaled, j * scaled, scaled);
      }
    }
  }
  
  getIndex(x, y) {
    //Check for valid x and y values
    x = max(min(x, this.n - 1), 0);
    y = max(min(y, this.n - 1), 0);
    
    return x + y * this.n;
  }
  
  addDensity(x, y, amount) {
    this.density[this.getIndex(x, y)] += amount;
  }
  
  addVelocity(x, y, amountX, amountY) {
    this.Vx[this.getIndex(x, y)] += amountX;
    this.Vy[this.getIndex(x, y)] += amountY;
  }
  
  setBound(b, x) {
    //Check ceiling and floor
    for(let i = 1; i < this.n - 1; i++) {
      x[this.getIndex(i, 0)] = b == 2 ? -x[this.getIndex(i, 1)] : x[this.getIndex(i, 1)];
      x[this.getIndex(i, this.n - 1)] = b == 2 ? -x[this.getIndex(i, this.n - 2)] : x[this.getIndex(i, this.n - 2)];
    }

    //Check walls
    for(let j = 1; j < this.n - 1; j++) {
      x[this.getIndex(0, j)] = b == 1 ? -x[this.getIndex(1, j)] : x[this.getIndex(1, j)];
      x[this.getIndex(this.n - 1, j)] = b == 1 ? -x[this.getIndex(this.n - 2, j)] : x[this.getIndex(this.n - 2, j)];
    }
    
    //Check corners
    x[this.getIndex(0, 0)] = 0.5 * (x[this.getIndex(1, 0)] + x[this.getIndex(0, 1)]);
    x[this.getIndex(0, this.n - 1)] = 0.5 * (x[this.getIndex(1, this.n - 1)] + x[this.getIndex(0, this.n - 2)]);
    x[this.getIndex(this.n - 1, 0)] = 0.5 * (x[this.getIndex(this.n - 2, 0)] + x[this.getIndex(this.n - 1, 1)]);
    x[this.getIndex(this.n - 1, this.n - 1)]   = 0.5 * (x[this.getIndex(this.n - 2, this.n - 1)] + x[this.getIndex(this.n - 1, this.n - 2)]);
  }

  solveLinear(b, x, x0, a, c) {
    //Use Gauss-Seidel method to approximate solution
    for (let k = 0; k < this.iter; k++) {
      for (let j = 1; j < this.n - 1; j++) {
        for (let i = 1; i < this.n - 1; i++) {
          x[this.getIndex(i, j)] = (x0[this.getIndex(i, j)] +
                                    a * (x[this.getIndex(i + 1, j)] +
                                    x[this.getIndex(i - 1, j)] +
                                    x[this.getIndex(i, j + 1)] +
                                    x[this.getIndex(i, j - 1)])) / c;
                                    
        }
      }
      this.setBound(b, x);
    }
  }
  
  diffuse(b, x, x0, diff) {
    //Apply diffusion
    let a = this.dt * diff * (this.n - 2) * (this.n - 2);
    this.solveLinear(b, x, x0, a, 1 + 6 * a);
  }
  
  project(Vx, Vy, p, div) {
    //Calculate net inflow/outflow
    for (let j = 1; j < this.n - 1; j++) {
      for (let i = 1; i < this.n - 1; i++) {
        div[this.getIndex(i, j)] = -0.5 * (Vx[this.getIndex(i + 1, j)] - 
                                           Vx[this.getIndex(i - 1, j)] + 
                                           Vy[this.getIndex(i, j + 1)] - 
                                           Vy[this.getIndex(i, j - 1)]) / this.n;
        p[this.getIndex(i, j)] = 0;
      }
    }
    
    //Calculate velocities required for net balance
    this.setBound(0, div); 
    this.setBound(0, p);
    this.solveLinear(0, p, div, 1, 6);
    
    for (let j = 1; j < this.n - 1; j++) {
      for (let i = 1; i < this.n - 1; i++) {
        Vx[this.getIndex(i, j)] -= 0.5 * (p[this.getIndex(i + 1, j)] - p[this.getIndex(i - 1, j)]) * this.n;
        Vy[this.getIndex(i, j)] -= 0.5 * (p[this.getIndex(i, j + 1)] - p[this.getIndex(i, j - 1)]) * this.n;
      }
    }
    
    this.setBound(1, Vx);
    this.setBound(2, Vy);
  }
  
  advect(b, d, d0, Vx, Vy) {
    let i0, i1, j0, j1;
    
    let dtx = this.dt * (this.n - 2);
    let dty = this.dt * (this.n - 2);
    
    let s0, s1, t0, t1;
    let tmp1, tmp2, x, y;
    
    //Predict previous values
    for(let j = 1; j < this.n - 1; j++) { 
      for(let i = 1; i < this.n - 1; i++) {
        tmp1 = dtx * Vx[this.getIndex(i, j)];
        tmp2 = dty * Vy[this.getIndex(i, j)];
        x = i - tmp1; 
        y = j - tmp2;
             
        x = max(x, 0.5);
        x = min(x, (this.n - 2) + 0.5);
        i0 = floor(x); 
        i1 = i0 + 1;
                
        y = max(y, 0.5);
        y = min(y, (this.n - 2) + 0.5);
        j0 = floor(y);
        j1 = j0 + 1; 
                
        s1 = x - i0; 
        s0 = 1 - s1; 
        t1 = y - j0; 
        t0 = 1 - t1;
                
        let i0i = parseInt(i0);
        let i1i = parseInt(i1);
        let j0i = parseInt(j0);
        let j1i = parseInt(j1);
        
        //Average values of surrounding cells
        d[this.getIndex(i, j)] = s0 * (t0 * d0[this.getIndex(i0i, j0i)] + t1 * d0[this.getIndex(i0i, j1i)]) +
                                 s1 * (t0 * d0[this.getIndex(i1i, j0i)] + t1 * d0[this.getIndex(i1i, j1i)]);
      }
    }
    
    this.setBound(b, d);
  }
}
