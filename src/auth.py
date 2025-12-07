from aiosqlite import Connection
from fastapi import Request, Response

import dotenv
import json
import jwt
import os

dotenv.load_dotenv()
JWT_SECRET = os.getenv("JWT_SECRET")
if JWT_SECRET is None or JWT_SECRET == "":
  raise ValueError("JWT_SECRET is not provided in the .env")

class User:
  name: str
  username: str

def gen_error(code: int, msg: str) -> Response:
  return Response(
    status_code=code,
    content=json.dumps({"error": msg}),
    headers={"Content-Type": "application/json"}
  )

def is_protected_path(path: str) -> bool:
  splits = path.split("/")
  return splits[1] == "discover" or splits[1] == "results"

async def validate_token(db: Connection, token: str | None) -> User | None:
  try:
    if token is None:
      return None

    decoded = jwt.decode(token, JWT_SECRET, algorithms="HS256")
    sub = decoded["sub"]
    if type(sub) is not str:
      return None
    if not sub.isnumeric():
      return None
    
    user = None
    cur = await db.execute("SELECT * FROM user WHERE id = ?", (int(sub),))
    rows = await cur.fetchall()
    for row in rows:
      user = User()
      user.name = row[1] # name field is 1st field
      user.username = row[2] # username field is 2nd field

    await cur.close()

    return user
  except Exception as e:
    print(e) # Debugging purpose
    return None

async def AuthMiddleware(request: Request, call_next):
  if not hasattr(request.app.state, "db"):
    print("App state doesn't have db attribute")
    return gen_error(500, "Internal server error")

  db = request.app.state.db
  if not isinstance(db, Connection):
    print("Invalid db state")
    return gen_error(500, "Internal server error")

  token = request.cookies.get("chem-token")
  validated = await validate_token(db, token)
  request.state.user = validated

  path = request.url.path

  if is_protected_path(path) and (validated is None):
    return gen_error(401, "Unauthorized")
  
  return await call_next(request)